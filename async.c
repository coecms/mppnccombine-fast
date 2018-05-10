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
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "netcdf.h"
#include "error.h"

#define TAG_CONTINUE       1
#define TAG_WRITE_CHUNK    2
#define TAG_OPEN_VARIABLE  3
#define TAG_CLOSE_VARIABLE 4
#define TAG_CLOSE          5

#define MAX_VARIABLES 100

typedef struct {
    hid_t file_id;
    
    struct {
        hid_t var_id;
        char varname[NC_MAX_NAME+1];
        size_t refcount;
    } vars[MAX_VARIABLES];
    
    int total_vars;
} async_state_t;

varid_t open_variable_async(
    const char * varname,
    size_t len,
    int async_writer_rank
    ) {
    varid_t out;

    MPI_Send(varname, len, MPI_CHAR, async_writer_rank,
             TAG_OPEN_VARIABLE, MPI_COMM_WORLD);
    MPI_Recv(&(out.idx), 1, MPI_INT, async_writer_rank,
             TAG_OPEN_VARIABLE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return out;
}

void close_variable_async(
    varid_t varid,
    int async_writer_rank
    ) {
    MPI_Send(&(varid.idx), 1, MPI_INT, async_writer_rank, TAG_CLOSE_VARIABLE, MPI_COMM_WORLD);
}

// Write a chunk in async mode
void write_chunk_async(
    varid_t var,
    size_t ndims,
    uint32_t filter_mask,
    hsize_t offset[],
    size_t data_size,
    const void * buffer,
    int async_writer_rank,
    MPI_Request * request
    ) {
    MPI_Request requests[4];

    MPI_Isend(&(var.idx), 1, MPI_INT, async_writer_rank,
              TAG_WRITE_CHUNK, MPI_COMM_WORLD, &(requests[0]));
    
    uint64_t ndims_ = ndims;
    MPI_Isend(&ndims_, 1, MPI_UINT64_T, async_writer_rank,
              TAG_CONTINUE, MPI_COMM_WORLD, &(requests[1]));
    MPI_Isend(&filter_mask, 1, MPI_UINT32_T, async_writer_rank,
              TAG_CONTINUE, MPI_COMM_WORLD, &(requests[2]));

    uint64_t offset_[ndims];
    for (int d=0; d<ndims; ++d) {offset_[d] = offset[d];}
    MPI_Isend(offset_, ndims, MPI_UINT64_T, async_writer_rank,
              TAG_CONTINUE, MPI_COMM_WORLD, &(requests[3]));

    MPI_Isend(buffer, data_size, MPI_CHAR, async_writer_rank,
              TAG_CONTINUE, MPI_COMM_WORLD, request);
}

static size_t receive_write_chunk_async(
    async_state_t * state,
    MPI_Status status
    ) {

    int idx;
    uint64_t ndims;
    uint32_t filter_mask;

    MPI_Recv(&idx, 1, MPI_INT, status.MPI_SOURCE,
             TAG_WRITE_CHUNK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ndims, 1, MPI_UINT64_T, status.MPI_SOURCE,
             TAG_CONTINUE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

    H5ERR(H5DOwrite_chunk(state->vars[idx].var_id, H5P_DEFAULT, filter_mask,
                          offset_, buffer_size, buffer));

    return buffer_size;
}


static void receive_open_variable_async(
    async_state_t * state,
    MPI_Status status
    ) {

    int len;
    MPI_Get_count(&status, MPI_CHAR, &len);

    char varname[len];
    MPI_Recv(varname, len, MPI_CHAR, status.MPI_SOURCE,
             TAG_OPEN_VARIABLE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int out = -1;
    for (int i=0; i<state->total_vars; ++i) {
        // See if the var exists
        if (strncmp(varname, state->vars[i].varname, len) == 0) {
            out = i;
            break;
        }
    }

    if (out < 0) {
        // Add a new var
        out = (state->total_vars)++;
        assert(state->total_vars < MAX_VARIABLES);
        assert(len <= NC_MAX_NAME + 1);

        state->vars[out].refcount = 0;
        strncpy(state->vars[out].varname, varname, len);
    }

    if (state->vars[out].refcount == 0) {
        // Open the var
        state->vars[out].var_id = 
            H5Dopen(state->file_id, state->vars[out].varname,
                    H5P_DEFAULT);
        H5ERR(state->vars[out].var_id);
    }

    // Increment the refcount
    state->vars[out].refcount++;

    // Send the id
    MPI_Send(&out, 1, MPI_INT, status.MPI_SOURCE, TAG_OPEN_VARIABLE,
             MPI_COMM_WORLD);

}

static void receive_close_variable_async(
    async_state_t * state,
    MPI_Status status) {

    int idx = 0;
    int comm_size;
    MPI_Recv(&idx, 1, MPI_INT, status.MPI_SOURCE, TAG_CLOSE_VARIABLE,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Decrement the reference count
    state->vars[idx].refcount--;

    // If zero has been reached close the variable
    if (state->vars[idx].refcount == 0) {
        H5ERR(H5Dclose(state->vars[idx].var_id));
    }
}

// Close the file (collective op)
void close_async(
    int async_writer_rank
    ) {

    MPI_Request r;
    MPI_Ibarrier(MPI_COMM_WORLD, &r);
    MPI_Wait(&r, MPI_STATUS_IGNORE);

    int buffer = 0;
    MPI_Send(&buffer, 1, MPI_INT, async_writer_rank, TAG_CLOSE, MPI_COMM_WORLD);
}

static void receive_close_async(
    async_state_t * state,
    MPI_Status status) {
    int buffer = 0;
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    for (int i=1;i<comm_size;++i) {
        MPI_Recv(&buffer, 1, MPI_INT, i, TAG_CLOSE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    H5ERR(H5Fclose(state->file_id));
}

// Async runner to accept writes
size_t run_async_writer(
    const char * filename
    ) {
    size_t total_size = 0;

    async_state_t state;

    state.file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    H5ERR(state.file_id);

    state.total_vars = 0;

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

            // Check for messages more often than barriers
            for (int j=0; j<10; ++j) {
                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
                if (flag) {
                    switch(status.MPI_TAG) {
                        case TAG_WRITE_CHUNK:
                            total_size += receive_write_chunk_async(&state, status);
                            break;
                        case TAG_OPEN_VARIABLE:
                            receive_open_variable_async(&state, status);
                            break;
                        case TAG_CLOSE_VARIABLE:
                            receive_close_variable_async(&state, status);
                            break;
                    }
                }
            }

            MPI_Test(&barrier_request, &reached_barrier, MPI_STATUS_IGNORE);
        }

        // Collective operations
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        switch(status.MPI_TAG) {
            case TAG_CLOSE:
                receive_close_async(&state, status);
                return total_size;
                break;
            default:
                fprintf(stderr, "Unknown TAG %d received by run_async_writer\n", status.MPI_TAG);
                MPI_Abort(MPI_COMM_WORLD, -1);
        }

    }
}
