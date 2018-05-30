/** 
 * Copyright 2018 
 *
 * \author   <scott.wales@unimelb.edu.au>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "error.h"

#include "hdf5.h"
#include "netcdf.h"
#include "mpi.h"

// NetCDF error handler
void handle_nc_error(int err, const char * file, int line) {
    if (err != 0) {
        const char * message = nc_strerror(err);

        fprintf(stderr, "ERROR in NetCDF %s:%d %d %s\n", file, line, err, message);
        MPI_Abort(MPI_COMM_WORLD, err);
    }
}

// HDF5 error handler
void handle_h5_error(int err, const char * file, int line) {
    if (err < 0) {
        fprintf(stderr, "ERROR in HDF5 %s:%d\n", file, line);
        H5Eprint1(stderr);
        MPI_Abort(MPI_COMM_WORLD, err);
    }
}

void handle_c_error(int err, const char * file, int line) {
    if (err != 0) {
        fprintf(stderr, "ERROR %s:%d\n", file, line);
        MPI_Abort(MPI_COMM_WORLD, err);
    }
}
