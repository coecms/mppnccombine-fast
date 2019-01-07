BUILD_DIR=build

all .DEFAULT:
	mkdir -p ${BUILD_DIR}
	cd ${BUILD_DIR} && cmake .. && ${MAKE} $@

clean:
	${RM} -r ${BUILD_DIR}

.SUFFIXES:
