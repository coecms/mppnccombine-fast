# mppnccombine-fast
[![DOI](https://zenodo.org/badge/131938571.svg)](https://zenodo.org/badge/latestdoi/131938571)

An accelerated version of the `mppnccombine` post-processing tool for MOM

Uses HDF5's raw IO functions to speed up collating large datasets - a 0.1
degree model goes from taking 4 hours to collate a compressed variable with
mppnccombine, to 6 minutes with mppnccombine-fast running with 16 processes

## Build

`mppnccombine-fast` requires HDF5 version 1.10.2 or above, as well as NetCDF 4,
a C compiler and a MPI library.

Cmake is used for building. A Makefile is also provided for ease of use:

    make # Release build

    make BUILD_TYPE=Debug # Debug build

    make check # Run tests

    make PREFIX=/apps/mppncc-fast install # Install mppnccombine-fast to $PREFIX

On Raijin the Makefile loads all required modules

The environment variables `$OPENMPI_ROOT`, `$HDF5_ROOT` and `$NETCDF_ROOT` may
be used to help locate libraries at other sites

## Use

Use like

    mpirun -n 2 ./mppnccombine-fast --output out.nc input.nc.0000 input.nc.0001 input.nc.0002

Files will be collated along all axes with a `domain_distribution` attribute

At least 2 MPI ranks need to be used (rank 0 writes the output file, other
ranks read). More can be used - input files will be balanced between the MPI
ranks.

## Commentary

The main slowdown in copying compressed variables is that the hdf5 library has
to de-compress them during the read, and re-compress them during the write.
`mppnccombine-fast` works around this by using HDF5 1.10.2's direct IO
functions
[H5DOwrite_chunk](https://support.hdfgroup.org/HDF5/doc/HL/RM_HDF5Optimized.html#H5DOwrite_chunk)
and
[H5DOread_chunk](https://support.hdfgroup.org/HDF5/doc/HL/RM_HDF5Optimized.html#H5DOread_chunk)
to copy the compressed data from one file to the other directly, rather than
going through the de-compress/re-compress cycle.

Since the NetCDF4 library is much nicer to use, but doesn't provide public
access to the underlying HDF5 file, we need to do a bit of musical chairs with
the files.

 1. The `init()` function
     1. Open the output file and the first input file in netcdf mode
     2. Copy NetCDF metadata and un-collated variables using the NetCDF library
     3. Close the NetCDF files
 2. The `copy()` function
     1. Open the output file in HDF5 mode
     2. For each input file:
         1. Open the input file in NetCDF mode
         2. Get the collated variables, sizes and offsets
         3. Re-open the input file in HDF5 mode
         4. Do a raw copy of the variables from the input to output files
         5. Close the input file
     3. Close the output file

To get a even larger speedup MPI is used to have separate read and write
processes, since HDF5 IO is a blocking function.

The communication between the read and write processes is handled by the file
`async.c` - the writer process runs a busy loop waiting for messages from the
reader processes, then handles messages as they come in. Individual reader
processes can be sending different variables at the same time.
