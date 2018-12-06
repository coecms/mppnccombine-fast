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
#include <stdbool.h>
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "netcdf.h"
#include "error.h"

//#define TAG_CONTINUE       1
#define TAG_CLOSE          5

#define TAG_WRITE_CHUNK    10
#define TAG_WRITE_CHUNK_NDIMS 11
#define TAG_WRITE_CHUNK_MASK 12
#define TAG_WRITE_CHUNK_OFFSET 13
#define TAG_WRITE_CHUNK_DATA 14

#define TAG_WRITE_FILTER   15
#define TAG_WRITE_FILTER_SHAPE   16
#define TAG_WRITE_FILTER_DATA   17

#define TAG_OPEN_VARIABLE  20
#define TAG_CLOSE_VARIABLE 21
#define TAG_VAR_INFO       22

#define MAX_VARIABLES 100
#define MAX_DIMS 10

typedef struct {
    hid_t file_id;
    
    struct {
        hid_t var_id;
        char varname[NC_MAX_NAME+1];
        size_t refcount;
        hsize_t chunk_shape[MAX_DIMS]; // Set by get_var_info
    } vars[MAX_VARIABLES];
    
    int total_vars;
} async_state_t;

varid_t open_variable_async(
    const char * varname,
    size_t len,
    int async_writer_rank
    ) {
    varid_t out;

    log_message(LOG_DEBUG, "SEND open variable %s", varname);
    MPI_Ssend(varname, len, MPI_CHAR, async_writer_rank,
             TAG_OPEN_VARIABLE, MPI_COMM_WORLD);
    log_message(LOG_DEBUG, "Waiting on %d %d", async_writer_rank, TAG_OPEN_VARIABLE);
    MPI_Recv(&(out.idx), 1, MPI_INT, async_writer_rank,
             TAG_OPEN_VARIABLE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    log_message(LOG_DEBUG, "RECV id of variable %s is %d", varname, out);

    return out;
}

// Get a dataspace, making sure it can fit the field
hid_t get_resized_space(async_state_t * state, int idx, int ndims, hsize_t offset[], hsize_t shape[]) {
    hid_t data_space = H5Dget_space(state->vars[idx].var_id);
    H5ERR(data_space);

    // Check if we need to extend the dataset
    hsize_t file_size[ndims];
    hsize_t max_size[ndims];
    H5ERR(H5Sget_simple_extent_dims(data_space, file_size, max_size));
    bool needs_resize = false;
    for (int d=0;d<ndims;++d) {
        if (offset[d] + shape[d] > file_size[d]) {
            needs_resize = true;
            file_size[d] = offset[d] + shape[d];
            if ((max_size[d] != H5S_UNLIMITED) && (file_size[d] > max_size[d])) {
                file_size[d] = max_size[d];
            }
        }
    }
    if (needs_resize) {
        H5ERR(H5Dset_extent(state->vars[idx].var_id, file_size));
        // Re-open
        data_space = H5Dget_space(state->vars[idx].var_id);
        H5ERR(data_space);
    }
    return data_space;
}

void close_variable_async(
    varid_t varid,
    int async_writer_rank
    ) {
    //log_message(LOG_DEBUG, "SEND close variable %d",varid);
    //MPI_Send(&(varid.idx), 1, MPI_INT, async_writer_rank, TAG_CLOSE_VARIABLE, MPI_COMM_WORLD);
}

// Change a NetCDF type to MPI
MPI_Datatype type_nc_to_mpi(nc_type type) {
    switch (type) {
        case (NC_INT):
            return MPI_INT;
        case (NC_FLOAT):
            return MPI_FLOAT;
        case (NC_DOUBLE):
            return MPI_DOUBLE;
        default:
            log_message(LOG_ERROR, "Unknown NetCDF type %d\n", type);
            MPI_Abort(MPI_COMM_WORLD, -1);
    }

    return MPI_INT;
}
// Change a NetCDF type to HDF5
hid_t type_nc_to_h5(nc_type type) {
    switch (type) {
        case (NC_INT):
            return H5T_NATIVE_INT;
        case (NC_FLOAT):
            return H5T_NATIVE_FLOAT;
        case (NC_DOUBLE):
            return H5T_NATIVE_DOUBLE;
        default:
            log_message(LOG_ERROR, "Unknown NetCDF type %d\n", type);
            MPI_Abort(MPI_COMM_WORLD, -1);
    }

    return H5T_NATIVE_INT;
}


// Get info about a variable
void variable_info_async(
    varid_t var,
    size_t ndims,
    uint64_t chunk[],
    int * deflate,
    int * deflate_level,
    int * shuffle,
    int async_writer_rank
    ) {

    // Send the variable ID we want to write to
    log_message(LOG_DEBUG, "SEND info variable %d", var);
    MPI_Ssend(&(var.idx), 1, MPI_INT, async_writer_rank,
             TAG_VAR_INFO, MPI_COMM_WORLD);

    uint64_t varinfo[ndims+3];
    MPI_Recv(varinfo, ndims+3, MPI_UINT64_T, async_writer_rank,
             TAG_VAR_INFO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    log_message(LOG_DEBUG, "RECV info variable %d", var);

    for (size_t d=0; d<ndims; ++d){
        chunk[d] = varinfo[d];
    }

    *deflate = varinfo[ndims];
    *deflate_level = varinfo[ndims+1];
    *shuffle = varinfo[ndims+2];
}

void receive_variable_info_async(
    async_state_t * state,
    MPI_Status status
    ) {

    int idx;
    MPI_Recv(&idx, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    log_message(LOG_DEBUG, "RECV info variable %s", state->vars[idx].varname);

    hid_t space = H5Dget_space(state->vars[idx].var_id);
    H5ERR(space);
    int ndims = H5Sget_simple_extent_ndims(space);
    H5ERR(H5Sclose(space));

    // Info is [ chunking[], deflate, deflate_level, shuffle ]
    unsigned long long varinfo[ndims + 3];
    hid_t plist = H5Dget_create_plist(state->vars[idx].var_id);
    H5ERR(plist);

    H5ERR(H5Pget_chunk(plist, ndims, varinfo + 0));
    H5ERR(H5Pget_chunk(plist, ndims, state->vars[idx].chunk_shape));

    int nfilters = H5Pget_nfilters(plist);
    H5ERR(nfilters);

    // Prep the filter info
    varinfo[ndims]   = 0;
    varinfo[ndims+1] = 0;
    varinfo[ndims+2] = 0;

    for (int filter_id=0; filter_id<nfilters; ++filter_id) {
        unsigned int flags;
        size_t elements = 4;
        unsigned values[elements];
        char name[256];
        unsigned config;
        H5Z_filter_t filter = H5Pget_filter(plist, filter_id, &flags,
                                            &elements, values, sizeof(name), name,
                                            &config);
        if (filter == H5Z_FILTER_DEFLATE) {
            // Store the deflate level
            varinfo[ndims] = 1;
            varinfo[ndims+1] = values[0];
        }
        if (filter == H5Z_FILTER_SHUFFLE) {
            // Store the shuffle
            varinfo[ndims+2] = 1;
        }

        H5ERR(filter);
    }

    H5ERR(H5Pclose(plist));

    MPI_Ssend(varinfo, ndims+3, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE,
             status.MPI_TAG, MPI_COMM_WORLD);
    log_message(LOG_DEBUG, "RECV info variable %s done", state->vars[idx].varname);
}


// Write data to the file, using the dataset filters
void write_uncompressed_async(
    varid_t var,
    size_t ndims,
    const size_t chunk_offset[],
    const size_t chunk_shape[],
    const void * buffer,
    nc_type type,
    int async_writer_rank,
    MPI_Request * request) {
    *request = MPI_REQUEST_NULL; 
    MPI_Request requests[4];

    // Send the variable ID we want to write to
    log_message(LOG_DEBUG, "SEND write variable %d", var);
    MPI_Ssend(&(var.idx), 1, MPI_INT, async_writer_rank,
              TAG_WRITE_FILTER, MPI_COMM_WORLD);

    // Pack chunk information into a vector
    int chunk_data_count = 2 * ndims + 1;
    uint64_t chunk_data[chunk_data_count];
    size_t buffer_count = 1;
    
    for (size_t d=0;d<ndims;++d) {
        chunk_data[d]       = chunk_offset[d];
        chunk_data[ndims+d] = chunk_shape[d];

        buffer_count *= chunk_shape[d];
    }

    chunk_data[2*ndims] = type;

    MPI_Ssend(chunk_data, chunk_data_count, MPI_UINT64_T, async_writer_rank,
              TAG_WRITE_FILTER_SHAPE, MPI_COMM_WORLD);


    MPI_Datatype type_mpi = type_nc_to_mpi(type);

    // Send the buffer - reciever can get the count and type
    MPI_Ssend(buffer, buffer_count, type_mpi, async_writer_rank,
              TAG_WRITE_FILTER_DATA, MPI_COMM_WORLD);

    int receipt;
    MPI_Recv(&receipt, 1, MPI_INT, async_writer_rank, TAG_WRITE_FILTER, 
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    log_message(LOG_DEBUG, "RECEIPT write variable %d", var);
}

static size_t receive_write_uncompressed_async(
    async_state_t * state,
    MPI_Status status
    ) {

    int idx;

    // Get the messages - target variable id, chunk info, data

    MPI_Recv(&idx, 1, MPI_INT, status.MPI_SOURCE,
             TAG_WRITE_FILTER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    log_message(LOG_DEBUG, "RECV write variable %s", state->vars[idx].varname);

    MPI_Status probe;
    MPI_Probe(status.MPI_SOURCE, TAG_WRITE_FILTER_SHAPE, MPI_COMM_WORLD, &probe);

    int chunk_data_count;
    MPI_Get_count(&probe, MPI_UINT64_T, &chunk_data_count);

    uint64_t chunk_data[chunk_data_count];
    MPI_Recv(chunk_data, chunk_data_count, MPI_UINT64_T, status.MPI_SOURCE,
             TAG_WRITE_FILTER_SHAPE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int type = chunk_data[chunk_data_count - 1];

    MPI_Datatype type_mpi = type_nc_to_mpi(type);

    MPI_Probe(status.MPI_SOURCE, TAG_WRITE_FILTER_DATA, MPI_COMM_WORLD, &probe);

    int buffer_count;
    MPI_Get_count(&probe, type_mpi, &buffer_count);

    char buffer[buffer_count * sizeof(double)];
    MPI_Recv(buffer, buffer_count, type_mpi, status.MPI_SOURCE,
             TAG_WRITE_FILTER_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Send receipt
    MPI_Send(&buffer_count, 1, MPI_INT, status.MPI_SOURCE, TAG_WRITE_FILTER, MPI_COMM_WORLD);
    log_message(LOG_DEBUG, "RECEIPT write variable %s", state->vars[idx].varname);

    // Unpack chunk info
    int ndims = (chunk_data_count - 1) / 2;
    hsize_t offset[ndims];
    hsize_t shape[ndims];
    for (int d=0;d<ndims;++d) {
        offset[d] = chunk_data[d];
        shape[d]  = chunk_data[ndims+d];
    }

    hsize_t block_size_start;
    herr_t err = H5Dget_chunk_storage_size(state->vars[idx].var_id, offset, &block_size_start);
    if (err < 0) block_size_start = 0;

    // Create selections and write out the data
    hid_t mem_space = H5Screate_simple(ndims, shape, NULL);
    H5ERR(mem_space);

    hid_t data_space = get_resized_space(state, idx, ndims, offset, shape);

    H5ERR(H5Sselect_hyperslab(data_space, H5S_SELECT_SET,
                              offset, NULL, shape, NULL));

    hid_t type_h5 = type_nc_to_h5(type);
    H5ERR(H5Dwrite(state->vars[idx].var_id, type_h5,
                   mem_space, data_space, H5P_DEFAULT, buffer));

    H5ERR(H5Sclose(mem_space));
    H5ERR(H5Sclose(data_space));

    hsize_t block_size_end;
    H5ERR(H5Dget_chunk_storage_size(state->vars[idx].var_id, offset, &block_size_end));

    // Return the change in the block size, in case of partial writes
    return block_size_end - block_size_start;
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
    *request = MPI_REQUEST_NULL;
    MPI_Request requests[4];

    log_message(LOG_DEBUG, "SEND raw write %d", var);
    MPI_Ssend(&(var.idx), 1, MPI_INT, async_writer_rank,
              TAG_WRITE_CHUNK, MPI_COMM_WORLD);
    
    uint64_t ndims_ = ndims;
    MPI_Ssend(&ndims_, 1, MPI_UINT64_T, async_writer_rank,
              TAG_WRITE_CHUNK_NDIMS, MPI_COMM_WORLD);
    MPI_Ssend(&filter_mask, 1, MPI_UINT32_T, async_writer_rank,
              TAG_WRITE_CHUNK_MASK, MPI_COMM_WORLD);

    uint64_t offset_[ndims];
    for (size_t d=0; d<ndims; ++d) {offset_[d] = offset[d];}
    MPI_Ssend(offset_, ndims, MPI_UINT64_T, async_writer_rank,
              TAG_WRITE_CHUNK_OFFSET, MPI_COMM_WORLD);

    log_message(LOG_DEBUG, "WAIT raw write %d", var);

    MPI_Ssend(buffer, data_size, MPI_CHAR, async_writer_rank,
              TAG_WRITE_CHUNK_DATA, MPI_COMM_WORLD);

    MPI_Wait(request, MPI_STATUS_IGNORE);
    log_message(LOG_DEBUG, "SEND raw write %d done", var);

    int receipt;
    MPI_Recv(&receipt, 1, MPI_INT, async_writer_rank, TAG_WRITE_CHUNK, 
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    log_message(LOG_DEBUG, "RECEIPT write variable %d", var);
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
             TAG_WRITE_CHUNK_NDIMS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&filter_mask, 1, MPI_UINT64_T, status.MPI_SOURCE,
             TAG_WRITE_CHUNK_MASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    uint64_t offset[ndims];
    MPI_Recv(offset, ndims, MPI_UINT64_T, status.MPI_SOURCE,
             TAG_WRITE_CHUNK_OFFSET, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    log_message(LOG_DEBUG, "RECV raw write variable %s", state->vars[idx].varname);

    int buffer_size;
    MPI_Status buffer_status;
    MPI_Probe(status.MPI_SOURCE, TAG_WRITE_CHUNK_DATA, MPI_COMM_WORLD, &buffer_status);
    MPI_Get_count(&buffer_status, MPI_CHAR, &buffer_size);

    char buffer[buffer_size];
    MPI_Recv(buffer, buffer_size, MPI_CHAR, status.MPI_SOURCE,
             TAG_WRITE_CHUNK_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    log_message(LOG_DEBUG, "RECV raw write variable %s passing to hdf5", state->vars[idx].varname);

    // Send receipt
    MPI_Send(&buffer_size, 1, MPI_INT, status.MPI_SOURCE, TAG_WRITE_CHUNK, MPI_COMM_WORLD);
    log_message(LOG_DEBUG, "RECEIPT write variable %s", state->vars[idx].varname);

    hsize_t offset_[ndims];
    for (size_t d=0; d<ndims; ++d) {offset_[d] = offset[d];}

    hid_t data_space = get_resized_space(state, idx, ndims, offset_, state->vars[idx].chunk_shape);
    H5ERR(H5Sclose(data_space));

    int err = (H5DOwrite_chunk(state->vars[idx].var_id, H5P_DEFAULT, filter_mask,
                          offset_, buffer_size, buffer));
    if (err < 0) {
        log_message(LOG_ERROR, "var %d %s from %d dims %zu [%zu,%zu,%zu]\n", idx, state->vars[idx].varname, status.MPI_SOURCE, ndims, offset[0], offset[1], offset[2]);
    }
    H5ERR(err);

    log_message(LOG_DEBUG, "RECV raw write variable %s done", state->vars[idx].varname);
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

    log_message(LOG_DEBUG, "RECV open variable %s", varname);

    int out = -1;
    for (int i=0; i<state->total_vars; ++i) {
        // See if the var exists
        if (strncmp(varname, state->vars[i].varname, len) == 0) {
            out = i;
            break;
        }
    }

    if (out < 0) {
        log_message(LOG_DEBUG, "Variable %s not in cache", varname);
        // Add a new var
        out = (state->total_vars)++;
        assert(state->total_vars < MAX_VARIABLES);
        assert(len <= NC_MAX_NAME + 1);

        state->vars[out].refcount = 0;
        strncpy(state->vars[out].varname, varname, len);
    }
    log_message(LOG_DEBUG, "Variable %s has id %d", varname, out);

    if (state->vars[out].refcount == 0) {
        // Open the var
        state->vars[out].var_id = 
            H5Dopen(state->file_id, state->vars[out].varname,
                    H5P_DEFAULT);
        H5ERR(state->vars[out].var_id);
    }

    // Increment the refcount
    state->vars[out].refcount++;
    log_message(LOG_DEBUG, "Variable %s has refcount %d", varname, state->vars[out].refcount);
    log_message(LOG_DEBUG, "Message %d %d %d", out, status.MPI_SOURCE, TAG_OPEN_VARIABLE);

    // Send the id
    MPI_Ssend(&out, 1, MPI_INT, status.MPI_SOURCE, TAG_OPEN_VARIABLE,
             MPI_COMM_WORLD);
    log_message(LOG_DEBUG, "OPEN sent");

}

static void receive_close_variable_async(
    async_state_t * state,
    MPI_Status status) {

    int idx = 0;
    MPI_Recv(&idx, 1, MPI_INT, status.MPI_SOURCE, TAG_CLOSE_VARIABLE,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    log_message(LOG_DEBUG, "RECV close variable %s", state->vars[idx].varname);

    // Decrement the reference count
    state->vars[idx].refcount--;

    //// If zero has been reached close the variable
    //if (state->vars[idx].refcount == 0) {
    //    H5ERR(H5Dclose(state->vars[idx].var_id));
    //}
}

// Close the file (collective op)
void close_async(
    int async_writer_rank
    ) {
    log_message(LOG_DEBUG, "SEND close file");
    MPI_Send(NULL, 0, MPI_INT, async_writer_rank, TAG_CLOSE, MPI_COMM_WORLD);
}

// Async runner to accept writes
size_t run_async_writer(
    const char * filename
    ) {
    MPI_Comm comm = MPI_COMM_WORLD;
    size_t total_size = 0;

    async_state_t state;

    state.file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    H5ERR(state.file_id);

    state.total_vars = 0;

    int open_senders;
    MPI_Comm_size(comm, &open_senders);
    --open_senders; // This process isn't sending

    while (open_senders > 0) {
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
        log_message(LOG_DEBUG, "PROBE tag %d from rank %d", status.MPI_TAG, status.MPI_SOURCE);

        switch(status.MPI_TAG) {
            case (TAG_OPEN_VARIABLE):
                receive_open_variable_async(&state, status);
                break;
            case (TAG_VAR_INFO):
                receive_variable_info_async(&state, status);
                break;
            case (TAG_WRITE_CHUNK):
                total_size += receive_write_chunk_async(&state, status);
                break;
            case (TAG_WRITE_FILTER):
                total_size += receive_write_uncompressed_async(&state, status);
                break;
            case (TAG_CLOSE_VARIABLE):
                receive_close_variable_async(&state, status);
                break;
            case (TAG_CLOSE):
                --open_senders;
                MPI_Recv(NULL, 0, MPI_INT, status.MPI_SOURCE, TAG_CLOSE, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                break;
        }
    }

    for (int v=0; v<state.total_vars; ++v) {
        H5ERR(H5Dclose(state.vars[v].var_id));
    }

    H5ERR(H5Fclose(state.file_id));
    log_message(LOG_DEBUG, "DONE run_async_writer");

    return total_size;
}
