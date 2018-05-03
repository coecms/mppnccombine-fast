CFLAGS=-std=c99 -Wall -check-pointers=rw -g -traceback -I/apps/hdf5/1.10.2/include

LDLIBS= -lhdf5 -lhdf5_hl

all: run_combo

run_h5collate: h5collate
	./h5collate
	touch $@

run_combo: combo
	./combo
	touch $@
