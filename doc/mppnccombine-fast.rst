API
===

async.h
-------

Contains the async Writer loop and functions that the Reader process can use to
communicate with it

Reader side
~~~~~~~~~~~

.. type:: varid_t

    Remote handle to a variable in the output file on the Writer process

.. doxygenfunction:: open_variable_async
.. doxygenfunction:: variable_info_async
.. doxygenfunction:: write_uncompressed_async
.. doxygenfunction:: write_chunk_async
.. doxygenfunction:: close_variable_async
.. doxygenfunction:: close_async

Writer side
~~~~~~~~~~~

.. doxygenfunction:: run_async_writer

error.h
-------

Functions for reporting errors from the various libraries used

.. doxygendefine:: NCERR
.. doxygendefine:: H5ERR
.. doxygendefine:: CERR

.. doxygenfunction:: set_log_level
.. doxygenfunction:: log_message

.. doxygendefine:: LOG_DEBUG
.. doxygendefine:: LOG_INFO
.. doxygendefine:: LOG_WARNING
.. doxygendefine:: LOG_ERROR

read_chunked.h
--------------

Functions the Readers use to read chunks from the input files and send them to
the Writer

.. doxygenfunction:: get_collated_dim_decomp
.. doxygenfunction:: get_collated_dim_len
.. doxygenfunction:: get_collation_info
.. doxygenfunction:: is_collated

.. doxygenfunction:: copy_chunked

