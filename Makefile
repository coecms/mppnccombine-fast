BUILD_DIR=build

ENV=
ifneq (,$(findstring raijin,${HOSTNAME}))
    ENV=module load cmake/3.12.2 intel-cc/2019.0.117 netcdf/4.6.1 hdf5/1.10.2 openmpi/3.1.3;
endif

all .DEFAULT:
	mkdir -p ${BUILD_DIR}
	${ENV} cd ${BUILD_DIR} && cmake ${CURDIR} && ${MAKE} $@

clean:
	${RM} -r ${BUILD_DIR}

.SUFFIXES:
