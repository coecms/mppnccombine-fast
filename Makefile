CC      = ${COMPILER_ENV} mpicc
#CFLAGS  = -std=c99 -Wall -Werror -check-pointers=rw -g -O0 -traceback
CFLAGS  = -std=c99 -Wall -Werror -g -O2
LDLIBS  = -lnetcdf -lhdf5_hl -lhdf5 -lm

ifneq ($(filter raijin%, ${HOSTNAME}),)
    # Setup the Raijin environment
    COMPILER_ENV = module purge; module load intel-cc openmpi/3.0.1 netcdf/4.6.1 hdf5/1.10.2;
    TEST_ENV     = module purge; module load conda openmpi/3.0.1;
endif

ifneq (${CONDA_BUILD},)
    # Load dependencies from conda
    # conda create -c conda-forge -n mppnccombine-fast-build gcc openmpi 'hdf5>=1.10.2' libnetcdf
    # conda create -c conda-forge -n mppnccombine-fast-test openmpi pytest xarray
    # Run 'conda activate' first so 'source' will work
    SHELL = bash
    COMPILER_ENV = source activate mppnccombine-fast-build;
    TEST_ENV = source activate mppnccombine-fast-test;
endif

ifneq (${PROFILE},)
    CFLAGS += -fprofile-arcs -ftest-coverage
endif

all: mppnccombine-fast

test: mppnccombine-fast
	${TEST_ENV} pytest test.py

mppnccombine-fast: async.o error.o readers.o

clean:
	${RM} mppnccombine-fast *.o
