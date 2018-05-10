CC      = mpicc
CFLAGS  = -std=c99 -Wall -Werror -check-pointers=rw -g -traceback
LDLIBS  = -lnetcdf -lhdf5_hl -lhdf5

# # Conda
# CFLAGS += -I${CONDA_PREFIX}/include
# LDFLAGS = -L${CONDA_PREFIX}/lib -Wl,-rpath=${CONDA_PREFIX}/lib

# Use `module load netcdf/4.6.1 hdf5/1.10.2`
with_module:
	module purge; module load intel-cc openmpi/2.1.1 netcdf/4.6.1 hdf5/1.10.2 scorep/3.1; ${MAKE} all

all: mppnccombine-fast

mppnccombine-fast: async.o error.o

clean:
	${RM} mppnccombine-fast *.o

scorep: with_module
	SCOREP_EXPERIMENT_DIRECTORY=scorep mpirun  -n 3 ./mppnccombine-fast --output /g/data/w35/saw562/test0.nc /short/v45/aek156/access-om2/archive/01deg_jra55_ryf/output243/ocean/ocean_temp_3hourly.nc.0000
