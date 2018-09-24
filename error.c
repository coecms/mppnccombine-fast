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
#define _GNU_SOURCE

#include "error.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "hdf5.h"
#include "netcdf.h"
#include "mpi.h"

#include <unistd.h>
#include <execinfo.h>

void print_backtrace(void) {
    void * stack[10];
    size_t size = backtrace(stack, 10);

    backtrace_symbols_fd(stack, size, STDERR_FILENO);
}

// NetCDF error handler
void handle_nc_error(int err, const char * file, int line) {
    if (err != 0) {
        const char * message = nc_strerror(err);

        log_message(LOG_ERROR, "ERROR in NetCDF %s:%d %d %s\n", file, line, err, message);
        print_backtrace();
        MPI_Abort(MPI_COMM_WORLD, err);
    }
}

// HDF5 error handler
void handle_h5_error(int err, const char * file, int line) {
    if (err < 0) {
        log_message(LOG_ERROR, "ERROR in HDF5 %s:%d\n", file, line);
        H5Eprint1(stderr);
        print_backtrace();
        MPI_Abort(MPI_COMM_WORLD, err);
    }
}

void handle_c_error(int err, const char * message, const char * file, int line) {
    if (err != 0) {
        log_message(LOG_ERROR, "ERROR %s:%d %s\n", file, line, message);
        print_backtrace();
        MPI_Abort(MPI_COMM_WORLD, err);
    }
}

static int log_level;

void set_log_level(int level) {
    log_level = level;
}

void log_message(int level, const char * message, ...) {
    if (level <= log_level) {
        int rank;
        va_list vargs;
        va_start(vargs, message);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        char * rendered;
        vasprintf(&rendered, message, vargs);
        fprintf(stderr, "[rank %03d] %s\n", rank, rendered);
        free(rendered);
        va_end(vargs);
    }
}
