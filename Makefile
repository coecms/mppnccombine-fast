CC      = mpicc
CFLAGS  = -std=c99 -Wall -Werror -check-pointers=rw -g -traceback
LDLIBS  = -lnetcdf -lhdf5_hl -lhdf5

# # Conda
# CFLAGS += -I${CONDA_PREFIX}/include
# LDFLAGS = -L${CONDA_PREFIX}/lib -Wl,-rpath=${CONDA_PREFIX}/lib

# Use `module load netcdf/4.6.1 hdf5/1.10.2`
with_module:
	module purge; module load intel-cc openmpi/1.10.7 netcdf/4.6.1 hdf5/1.10.2; ${MAKE} test

test: mppnccombine-fast
	time mpiexec -n 2 ./mppnccombine-fast --output /g/data/${PROJECT}/${USER}/test.nc /short/v45/aek156/access-om2/archive/01deg_jra55_ryf/output243/ocean/ocean_temp_3hourly.nc.0000
	module purge; module load conda; ./view.py /g/data/${PROJECT}/${USER}/test.nc

all: mppnccombine-fast
