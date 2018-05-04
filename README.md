# mppnccombine-fast

An accelerated version of the `mppnccombine` post-processing tool for MOM

Uses HDF5's raw IO functions to speed up collating large datasets - a 0.1
degree model goes from taking 4 hours to collate a compressed variable with
mppnccombine, to 25 minutes with mppnccombine-fast

## Build

`mppnccombine-fast` requires HDF5 version 1.10.2 or above

On Raijin:

    module load netcdf/4.6.1 hdf5/1.10.2 intel-cc openmpi
    make

## Use

Use like

    mppnccombine-fast --output out.nc input.nc.0000 input.nc.0001 input.nc.0002

Files will be collated along all axes with a `domain_distribution` attribute


TODO: MPI can also be used

    mpirun mppnccombine-fast --output out.nc input.nc.0000 input.nc.0001 input.nc.0002

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
