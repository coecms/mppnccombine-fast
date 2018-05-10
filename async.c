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

#include "async.h"

#include <stdint.h>
#include <stdlib.h>
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "error.h"

#define TAG_CONTINUE       1
#define TAG_WRITE_CHUNK    2
#define TAG_OPEN_VARIABLE  3
#define TAG_CLOSE_VARIABLE 4
#define TAG_CLOSE          5

// Write a chunk in async mode to the currently open variable
void write_chunk_async(
    size_t ndims,
    uint32_t filter_mask,
    hsize_t offset[],
    size_t data_size,
    const void * buffer,
    int async_writer_rank,
    MPI_Request * request
    ) {
    MPI_Request requests[3];
    
    uint64_t ndims_ = ndims;
    MPI_Isend(&ndims_, 1, MPI_UINT64_T, async_writer_rank,
              TAG_WRITE_CHUNK, MPI_COMM_WORLD, &(requests[0]));
    MPI_Isend(&filter_mask, 1, MPI_UINT32_T, async_writer_rank,
              TAG_CONTINUE, MPI_COMM_WORLD, &(requests[1]));

    uint64_t offset_[ndims];
    for (int d=0; d<ndims; ++d) {offset_[d] = offset[d];}
    MPI_Isend(offset_, ndims, MPI_UINT64_T, async_writer_rank,
              TAG_CONTINUE, MPI_COMM_WORLD, &(requests[2]));

    MPI_Isend(buffer, data_size, MPI_CHAR, async_writer_rank,
              TAG_CONTINUE, MPI_COMM_WORLD, request);
}

void receive_write_chunk_async(
    hid_t var,
    MPI_Status status
    ) {

    uint64_t ndims;
    uint32_t filter_mask;

    MPI_Recv(&ndims, 1, MPI_UINT64_T, status.MPI_SOURCE,
             TAG_WRITE_CHUNK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&filter_mask, 1, MPI_UINT64_T, status.MPI_SOURCE,
             TAG_CONTINUE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    uint64_t offset[ndims];
    MPI_Recv(offset, ndims, MPI_UINT64_T, status.MPI_SOURCE,
             TAG_CONTINUE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int buffer_size;
    MPI_Status buffer_status;
    MPI_Probe(status.MPI_SOURCE, TAG_CONTINUE, MPI_COMM_WORLD, &buffer_status);
    MPI_Get_count(&buffer_status, MPI_CHAR, &buffer_size);

    char buffer[buffer_size];
    MPI_Recv(buffer, buffer_size, MPI_CHAR, status.MPI_SOURCE,
             TAG_CONTINUE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    hsize_t offset_[ndims];
    for (int d=0; d<ndims; ++d) {offset_[d] = offset[d];}

    H5ERR(H5DOwrite_chunk(var, H5P_DEFAULT, filter_mask, offset_, buffer_size, buffer));
}


// Open a variable in async mode (collective)
void open_variable_async(
    const char * varname,
    size_t len, // Total length including '/0'
    int async_writer_rank
    ) {

    MPI_Request r;
    MPI_Ibarrier(MPI_COMM_WORLD, &r);
    MPI_Wait(&r, MPI_STATUS_IGNORE);
    printf("Passed barrier\n");
    MPI_Send(varname, len, MPI_CHAR, async_writer_rank, TAG_OPEN_VARIABLE, MPI_COMM_WORLD);
}

void receive_open_variable_async(
    hid_t file_in,
    hid_t * var_out,
    MPI_Status status
    ) {

    int len;
    MPI_Get_count(&status, MPI_CHAR, &len);

    char varname[len];

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    for (int i=1;i<comm_size;++i) {
        MPI_Recv(varname, len, MPI_CHAR, i, TAG_OPEN_VARIABLE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    *var_out = H5Dopen(file_in, varname, H5P_DEFAULT);
}

void close_variable_async(
    int async_writer_rank
    ) {

    MPI_Request r;
    MPI_Ibarrier(MPI_COMM_WORLD, &r);
    MPI_Wait(&r, MPI_STATUS_IGNORE);

    int buffer = 0;
    MPI_Send(&buffer, 1, MPI_INT, async_writer_rank, TAG_CLOSE_VARIABLE, MPI_COMM_WORLD);

}

void receive_close_variable_async(hid_t var_in, MPI_Status status) {
    int buffer = 0;
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    for (int i=1;i<comm_size;++i) {
        MPI_Recv(&buffer, 1, MPI_INT, i, TAG_CLOSE_VARIABLE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    H5ERR(H5Dclose(var_in));
}

void close_async(
    int async_writer_rank
    ) {

    MPI_Request r;
    MPI_Ibarrier(MPI_COMM_WORLD, &r);
    MPI_Wait(&r, MPI_STATUS_IGNORE);

    int buffer = 0;
    MPI_Send(&buffer, 1, MPI_INT, async_writer_rank, TAG_CLOSE, MPI_COMM_WORLD);
}

void receive_close_async(hid_t file_in, MPI_Status status) {
    int buffer = 0;
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    for (int i=1;i<comm_size;++i) {
        MPI_Recv(&buffer, 1, MPI_INT, i, TAG_CLOSE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    H5ERR(H5Fclose(file_in));
}

// Async runner to accept writes
void run_async_writer(
    const char * filename
    ) {

    hid_t file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    H5ERR(file);

    hid_t var;

    while (true) {
        // Check for a collective operation
        MPI_Request barrier_request;
        MPI_Ibarrier(MPI_COMM_WORLD, &barrier_request);

        int reached_barrier = 0;
        MPI_Test(&barrier_request, &reached_barrier, MPI_STATUS_IGNORE);

        // Non-collective operations
        while (!reached_barrier) {
            MPI_Status status;
            int flag;

            // TAG_WRITE_CHUNK
            MPI_Iprobe(MPI_ANY_SOURCE, TAG_WRITE_CHUNK, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                receive_write_chunk_async(var, status);
            }

            MPI_Test(&barrier_request, &reached_barrier, MPI_STATUS_IGNORE);
        }

        // Collective operations
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        switch(status.MPI_TAG) {
            case TAG_OPEN_VARIABLE:
                receive_open_variable_async(file, &var, status);
                break;
            case TAG_CLOSE_VARIABLE:
                receive_close_variable_async(var, status);
                break;
            case TAG_CLOSE:
                receive_close_async(file, status);
                return;
                break;
            default:
                fprintf(stderr, "Unknown TAG %d received by run_async_writer\n", status.MPI_TAG);
                MPI_Abort(MPI_COMM_WORLD, -1);
        }

    }
}
