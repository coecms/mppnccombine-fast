CC      = module purge; module load intel-cc openmpi/3.0.1 netcdf/4.6.1 hdf5/1.10.2; mpicc
#CFLAGS  = -std=c99 -Wall -Werror -check-pointers=rw -g -traceback
CFLAGS  = -std=c99 -Wall -Werror -g -O2
LDLIBS  = -lnetcdf -lhdf5_hl -lhdf5

# # Conda
# CFLAGS += -I${CONDA_PREFIX}/include
# LDFLAGS = -L${CONDA_PREFIX}/lib -Wl,-rpath=${CONDA_PREFIX}/lib

all: mppnccombine-fast

test: mppnccombine-fast
	module purge; module load conda openmpi/3.0.1; pytest test.py


mppnccombine-fast: async.o error.o

clean:
	${RM} mppnccombine-fast *.o

scorep: CC_ = scorep mpicc
scorep: with_module
	SCOREP_EXPERIMENT_DIRECTORY=scorep mpirun -n 2 ./mppnccombine-fast --output /g/data/w35/saw562/test0.nc /short/v45/aek156/access-om2/archive/01deg_jra55_ryf/output243/ocean/ocean_temp_3hourly.nc.0000
