Implementation Overview
=======================

The basic outline of ``mppnccombine-fast`` consists of one "Writer" rank and one or
more "Reader" ranks. The Writer rank handles all writing to the output file,
while the Reader ranks read in data from the many files to be collated and
send the data to the Writer rank.

The main slowdown in copying compressed variables is that the hdf5 library has
to de-compress them during the read, and re-compress them during the write.
``mppnccombine-fast`` works around this by using HDF5 1.10.2's direct IO
functions
`H5DOwrite_chunk() <https://support.hdfgroup.org/HDF5/doc/HL/RM_HDF5Optimized.html#H5DOwrite_chunk>`_
and
`H5DOread_chunk() <https://support.hdfgroup.org/HDF5/doc/HL/RM_HDF5Optimized.html#H5DOread_chunk>`_
to copy the compressed data from one file to the other directly, rather than
going through the de-compress/re-compress cycle.

To get a even larger speedup MPI is used to have separate read and write
processes, since HDF5 IO is a blocking function.

Since the NetCDF4 library is much nicer to use, but doesn't provide public
access to the underlying HDF5 file, we need to do a bit of musical chairs with
the files, swapping between NetCDF4 and HDF5 modes by re-opening the files.

.. graphviz::
    
    digraph ranks {
        subgraph cluster_writer {
            copy_uncollated -> async_loop -> all_ended

            label = "Writer"
            copy_uncollated [label = "Copy Uncollated Vars"]
            async_loop [label = "Async Write Loop"]
            all_ended [label = "All Readers Ended"]
            rank = same
        }

        subgraph cluster_reader {
            open_file -> send_chunks -> end_of_files

            label = "Reader"
            open_file [label = "Open File"]
            send_chunks [label = "Send Collated Chunks"]
            end_of_files [label = "All Files Read"]
            rank = same
        }

        send_chunks:s -> open_file:e [constraint=false]

        async_loop:s -> async_loop:n [constraint=false]

        send_chunks -> async_loop [constraint=false]
        end_of_files -> async_loop [constraint=false]
    }

Writer Rank
-----------

The Writer starts out by copying the dimensions, attributes and any uncollated
variables from the first of the listed input files using the NetCDF API in
:func:`init()` and :func:`copy_contiguous()`. It then re-opens the file using
the HDF5 API and enters the 'Async Write Loop' in :func:`run_async_writer()`.

This loop polls for any incoming MPI messages from the Reader processes then
performs some action (e.g. write a compressed chunk directly to the file at
some location). Once a Reader has finished reading all of its input files it
sends a close message to the Writer rank, once all close messages have been
received the Writer rank closes the output file and exits.

Reader Ranks
------------

The Readers distribute input files amonst themselves using a shared atomic
counter. When a Reader is ready for a new file it gets the next value from the
counter, then in :func:`copy_chunked()` opens that file using NetCDF to query
its attributes and discover and copy collated variables.

Depending on the chunking and alignment of the file the Reader will decide to
copy the data either in uncompressed form using NetCDF with
:func:`copy_netcdf_variable_chunks()` or directly copying the compressed chunks
by re-opening the file in HDF5 mode with :func:`copy_hdf5_variable_chunks()`.

Once all available files have been read the Reader sends a final close message
to the Writer and exits.

The Async Write Loop
--------------------

The Async write loop is set up to handle a number of messages that the Readers
will send to the Writer

 * :cpp:func:`open_variable_async`: Obtain a handle to a variable in the output
   file
 * :func:`variable_info_async()`: Obtain chunking and compression information
   for a variable
 * :func:`write_uncompressed_async()`: Write uncompressed data to a given
   logical location in the variable
 * :func:`write_chunk_async()`: Write a compressed chunk directly to the output
   file at a given chunk location
 * :func:`close_variable_async()`: Return the variable handle
 * :func:`close_async()`: Reports that the Reader will not send any more
   messages

The Writer asyncronously polls for these messages in
:func:`run_async_writer()`, then actions them in
:func:`receive_open_variable_async()` etc.
