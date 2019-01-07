cmake_minimum_required(VERSION 3.10)
project(mppnccombine-fast)
include(FindMPI)

add_executable(mppnccombine-fast mppnccombine-fast.c async.c error.c read_chunked.c)


find_path(HDF5_INCLUDE_DIRECTORY
	NAMES hdf5.h
	HINTS ENV HDF5_ROOT
	PATH SUFFIXES include
	)
find_library(HDF5_LIB
	NAMES hdf5
	HINTS ENV HDF5_ROOT
	PATH SUFFIXES lib
	)
find_library(HDF5_HL_LIB
	NAMES hdf5_hl
	HINTS ENV HDF5_ROOT
	PATH SUFFIXES lib
	)
find_path(NETCDF_INCLUDE_DIRECTORY
	NAMES netcdf.h
	HINTS ENV NETCDF_ROOT
	PATH SUFFIXES include
	)
find_library(NETCDF_LIB
	NAMES netcdf
	HINTS ENV NETCDF_ROOT
	PATH SUFFIXES lib
	)

find_program(CLANG_TIDY
	NAMES clang-tidy
	)
if (CLANG_TIDY)
	set_target_properties(mppnccombine-fast PROPERTIES
		C_CLANG_TIDY "${CLANG_TIDY};-checks=*")
endif()

target_compile_definitions(mppnccombine-fast PRIVATE ${MPI_C_COMPILE_DEFINITIONS})
target_compile_options(mppnccombine-fast PRIVATE ${MPI_C_COMPILE_OPTIONS})
target_link_libraries(mppnccombine-fast PUBLIC
	-lm
	${MPI_C_LIBRARIES}
	${HDF5_LIB}
	${HDF5_HL_LIB}
	${NETCDF_LIB})
target_include_directories(mppnccombine-fast PRIVATE
	${MPI_C_HEADER_DIR}
	${HDF5_INCLUDE_DIRECTORY}
	${NETCDF_INCLUDE_DIRECTORY})
