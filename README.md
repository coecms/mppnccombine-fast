# mppnccombine-fast

An accelerated version of the `mppnccombine` post-processing tool for MOM

Uses HDF5's raw IO functions to speed up collating large datasets

## Build

`mppnccombine-fast` requires HDF5 version 1.10.2 or above

On Raijin:

    module load netcdf/4.6.1 hdf5/1.10.2 intel-cc openmpi
    make

## Use

TODO: Use like

    mppnccombine-fast out.nc input.nc.0000 input.nc.0001 input.nc.0002

Files will be collated along the `xt_ocean` and `yt_ocean` dimensions


TODO: MPI can also be used

   mpirun mppnccombine-fast out.nc input.nc.0000 input.nc.0001 input.nc.0002

