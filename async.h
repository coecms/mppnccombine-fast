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

#ifndef ASYNC_H
#define ASYNC_H
#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include "hdf5.h"
#include "mpi.h"

// Write a chunk in async mode to the currently open variable
void write_chunk_async(
    size_t ndims,
    uint32_t filter_mask,
    hsize_t offset[],
    size_t data_size,
    const void * buffer,
    int async_writer_rank,
    MPI_Request * request
    );

// Open a variable in the async writer (collective)
void open_variable_async(
    const char * varname,
    size_t len, // Total length including '/0'
    int async_writer_rank
    );

// Close a variable in the async writer (collective)
void close_variable_async(
    int async_writer_rank
    );

// Close the async writer (collective)
void close_async(
    int async_writer_rank
    );

// Async runner to accept writes
void run_async_writer(
    const char * filename
    );

#ifdef __cplusplus
}
#endif
#endif
