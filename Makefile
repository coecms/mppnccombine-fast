CFLAGS  = -std=c99 -Wall -check-pointers=rw -g -traceback
LDLIBS  = -lnetcdf -lhdf5_hl -lhdf5

# # Conda
# CFLAGS += -I${CONDA_PREFIX}/include
# LDFLAGS = -L${CONDA_PREFIX}/lib -Wl,-rpath=${CONDA_PREFIX}/lib
#
# Use `module load netcdf/4.6.1 hdf5/1.10.2`
#

all: run_combo

run_h5collate: h5collate
	./h5collate
	touch $@

run_combo: combo
	./combo
	touch $@
