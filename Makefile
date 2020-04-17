BUILD_DIR=build
BUILD_TYPE?=Release

CMAKE_FLAGS=-DCMAKE_BUILD_TYPE=${BUILD_TYPE}

ifdef PREFIX
    CMAKE_FLAGS+= -DCMAKE_INSTALL_PREFIX=${PREFIX}
endif

ifneq (,$(findstring gadi,${HOSTNAME}))
    # raijin.nci.org.au build environment
    BUILD_ENV=module purge; module load cmake/3.16.2 intel-compiler/2020.0.166 netcdf/4.7.3 hdf5/1.10.5 openmpi/4.0.2;
    TEST_ENV=module purge; module load conda openmpi/4.0.2;
endif

all .DEFAULT:
	mkdir -p ${BUILD_DIR}
	${BUILD_ENV} cd ${BUILD_DIR} && cmake ${CMAKE_FLAGS} ${CURDIR} && ${MAKE} $@

clean:
	${RM} -r ${BUILD_DIR}

check test: mppnccombine-fast
	${TEST_ENV} PATH=${BUILD_DIR}:$$PATH py.test

doc:
	doxygen
	${MAKE} -C doc html

.SUFFIXES:
.PHONY: all clean check test doc
