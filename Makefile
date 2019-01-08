BUILD_DIR=build
BUILD_TYPE?=Release

CMAKE_FLAGS=-DCMAKE_BUILD_TYPE=${BUILD_TYPE}

ifdef PREFIX
    CMAKE_FLAGS+= -DCMAKE_INSTALL_PREFIX=${PREFIX}
endif

ifneq (,$(findstring raijin,${HOSTNAME}))
    # raijin.nci.org.au build environment
    BUILD_ENV=module purge; module load cmake/3.12.2 intel-cc/2019.0.117 netcdf/4.6.1 hdf5/1.10.2 openmpi/3.0.1;
    TEST_ENV=module purge; module load conda openmpi/3.0.1;
endif

all .DEFAULT:
	mkdir -p ${BUILD_DIR}
	${BUILD_ENV} cd ${BUILD_DIR} && cmake ${CMAKE_FLAGS} ${CURDIR} && ${MAKE} $@

clean:
	${RM} -r ${BUILD_DIR}

check test: mppnccombine-fast
	${TEST_ENV} PATH=${BUILD_DIR}:$$PATH py.test

.SUFFIXES:
