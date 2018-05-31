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

#include "readers.h"
#include "error.h"
#include "async.h"

#include "netcdf.h"
#include "hdf5.h"
#include "hdf5_hl.h"

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

// Get the decomposition attribute from a variable
// If this is not a decomposed variable, decompositon[] is unchanged and returns false,
// otherwise returns true
bool get_collated_dim_decomp(int ncid, const char * varname, int decomposition[4]) {
    int varid;
    int err;
    NCERR(nc_inq_varid(ncid, varname, &varid));

    err = nc_get_att_int(ncid, varid, "domain_decomposition", decomposition);
    if (err == NC_ENOTATT) {
        return false;
    }
    NCERR(err);
    return true;
}

// Get the total length of a collated variable
void get_collated_dim_len(int ncid, const char * varname, size_t * len) {
    int decomposition[4];
    get_collated_dim_decomp(ncid, varname, decomposition);
    *len = decomposition[1];
}

// Get collation info from a variable
// in_offset[ndims]  - The offset in the local data TODO: remove
// out_offset[ndims] - The offset in the collated array of this variable
// local_size[ndims] - The size of this variable
// total_size[ndims] - The total collated size of this variable
// returns true if any of the dimensions are collated
bool get_collation_info(int ncid, int varid,
                        size_t in_offset[], size_t out_offset[],
                        size_t local_size[], size_t total_size[],
                        int ndims) {

    // Get the dimension ids
    int dimids[ndims];
    NCERR(nc_inq_vardimid(ncid, varid, dimids));

    bool is_collated[ndims];
    bool out = false;

    for (int d=0; d<ndims; ++d) {
        // Dimension name
        char dimname[NC_MAX_NAME+1];
        NCERR(nc_inq_dimname(ncid, dimids[d], dimname));

        // Get the decomposition
        int decomposition[4];
        is_collated[d] = get_collated_dim_decomp(ncid, dimname, decomposition);

        // Calculate the per-dim values
        if (is_collated[d]) {
            in_offset[d] = 0;
            out_offset[d] = decomposition[2] - 1;
            local_size[d] = decomposition[3] - out_offset[d];
            total_size[d] = decomposition[1];
        } else {
            in_offset[d] = 0;
            out_offset[d] = 0;
            size_t len;
            NCERR(nc_inq_dimlen(ncid, dimids[d], &len));
            local_size[d] = len;
            total_size[d] = len;
        }

        // Will be true if any dimension is collated
        out = out || is_collated[d];
    }
    return out;
}

// Returns true if any of the dimensions are collated
bool is_collated(int ncid, int varid) {
    int ndims;
    NCERR(nc_inq_varndims(ncid, varid, &ndims));

    // We'll discard these arrays
    size_t in_offset[ndims];
    size_t out_offset[ndims];
    size_t local_size[ndims];
    size_t total_size[ndims];

    return get_collation_info(ncid, varid, in_offset, out_offset, local_size, total_size, ndims);
}


struct chunker_t {
    hid_t file;
    hid_t dataset;
    int ndims;
    size_t chunk[10];
    void * buffer;
    size_t buffer_size;
};

struct chunker_props {
    int data_type;
    int ndims;
    size_t chunk[10];
    bool shuffle;
    int deflate_level;
    bool scale_offset;
    int scale_type;
    int scale_factor;
};

// Create a temporary dataset to handle compression
void chunker_create(struct chunker_t * chunker, struct chunker_props * props) {
    // We'll open this path exclusively using HDF5 - in the case of a race HDF5 will error
    char * temp_filename = "x";// = tempnam(NULL, "mnccf");

    // Open the file
    chunker->file = H5Fcreate(temp_filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    H5ERR(chunker->file);

    // Drop the file so it's deleted at the end
    CERR(unlink(temp_filename));

    // Make a dataspace
    hsize_t hchunk[props->ndims];
    for (int i=0;i<props->ndims;++i) hchunk[i] = props->chunk[i];
    hid_t space = H5Screate_simple(props->ndims, hchunk, NULL);
    H5ERR(space);

    // Make a dataset in the file
    hid_t dcpl = H5Pcreate(H5P_DATATYPE_CREATE);
    H5ERR(dcpl);
    H5ERR(H5Pset_chunk(dcpl, props->ndims, hchunk));

    chunker->dataset = H5Dcreate(chunker->file, "chunker", props->data_type, space,
                                 H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5ERR(chunker->dataset);

    H5ERR(H5Sclose(space));
    H5ERR(H5Pclose(dcpl));
}

// Close the temporary dataset
void chunker_close(struct chunker_t * chunker) {
    H5ERR(H5Dclose(chunker->dataset));
    H5ERR(H5Fclose(chunker->file));
    free(chunker->buffer);
}

// Read a hdf5 processed chunk, putting the data of size `data_size` into `buffer`
// If needed the buffer will be resized
void read_hdf5_chunk(hid_t dataset,
                     hsize_t offset[],
                     void ** buffer,
                     size_t * buffer_size,
                     size_t * data_size) {
    uint32_t mask;
    hsize_t compressed_size;
    H5ERR(H5Dget_chunk_storage_size(dataset, offset, &compressed_size));
    if (compressed_size > *buffer_size) {
        *buffer_size = compressed_size;
        *buffer = realloc(*buffer, *buffer_size);
    }
    H5ERR(H5DOread_chunk(dataset, H5P_DEFAULT, offset, &mask, *buffer));
    *data_size = compressed_size;
}

// Given information about the local data and chunking, get the offset and
// shape of a specific chunk, taking into account partial chunks that don't
// line up with the output file chunking
void get_chunk_offset_shape(size_t ndims,                // Number of dims
                            const size_t local_offset[], // Local data offset
                            const size_t local_shape[],  // Local data shape
                            const size_t chunk[],        // Chunk shape
                            const size_t nchunks[],      // Number of chunks in each dim
                            size_t c,                    // Chunk id
                            size_t chunk_offset_in[],    // Offset of the chunk local data
                            size_t chunk_offset_out[],   // Offset of the chunk remote data
                            size_t chunk_shape[],        // Shape of the chunk local data
                            bool * partial               // True if partial chunk
                           ) {

    *partial = false;

    for (int i=ndims-1; i>=0; --i) {
        chunk_offset_in[i] = -1;
        chunk_offset_out[i] = -1;
        chunk_shape[i] = -1;

        size_t chunkid = c % nchunks[i];

        // Start with a full sized chunk
        chunk_shape[i] = chunk[i];

        // Count from the chunk boundary to the left of the starting point
        chunk_offset_out[i] = floor(local_offset[i] / chunk[i])*chunk[i]
            + chunk[i] * chunkid;

        if (chunk_offset_out[i] < local_offset[i]) {
            // Partial chunk on the left
            *partial = true;
            chunk_shape[i] = chunk[i] - (local_offset[i] - chunk_offset_out[i]);
            chunk_offset_out[i] = local_offset[i];
        }
        printf("%zu %zu %zu %zu\n", c, chunk_offset_in[0], chunk_offset_out[0], chunk_shape[0]);

        // Total size in dimension 'i'
        ssize_t total_size = (chunk[i] - local_offset[i]%chunk[i])
            + chunkid * chunk[i];

        if (total_size > local_shape[i]) {
            // Partial chunk on the right
            *partial = true;
            chunk_shape[i] -= total_size - local_shape[i];
        }

        chunk_offset_in[i] = chunk_offset_out[i] - local_offset[i];

        c = c / nchunks[i];

        // Check we're in bounds of the input data
        assert(chunk_offset_in[i] < local_shape[i]);
        assert(chunk_offset_in[i] + chunk_shape[i] <= local_shape[i]);
    }
}


// Read `chunk` sized chunks from variable `varid` of netcdf file `ncid`,
// sending them to variable `var` of the async writer at `async_writer_rank`
void copy_netcdf_variable_chunks(
                                 varid_t var,
                                 int ncid,
                                 int varid,
                                 int ndims,
                                 const size_t chunk[],
                                 int async_writer_rank
                                ) {
    size_t local_offset[ndims];
    size_t shape[ndims];
    size_t nchunks[ndims];
    size_t total_chunks = 1;
    size_t chunk_size = 1;

    {
        size_t in_offset[ndims];
        size_t total_shape[ndims];

        get_collation_info(ncid, varid, in_offset, local_offset,
                           shape, total_shape, ndims);
    }

    int dimids[ndims];
    NCERR(nc_inq_vardimid(ncid, varid, dimids));

    for (int i=0;i<ndims;++i) {
        NCERR(nc_inq_dimlen(ncid, dimids[i], &(shape[i])));

        // Need to account for partial chunks at the start of this file's block
        nchunks[i] = ceil(((local_offset[i] % chunk[i]) + shape[i]) / (float)chunk[i]);
        printf("chunks %zu\n",nchunks[i]);

        total_chunks *= nchunks[i];
        chunk_size *= chunk[i];
    }

    void * buffer = malloc(chunk_size*sizeof(double));
    nc_type type;
    NCERR(nc_inq_vartype(ncid, varid, &type));

    for (size_t c=0; c<total_chunks; ++c) {
        size_t chunk_offset_in[ndims]; // Offset in the data to read from
        size_t chunk_offset_out[ndims]; // Offset in the data to write to
        size_t chunk_shape[ndims];  // Shape of this chunk
        bool partial_chunk;         // Is this chunk partial? (so send uncompressed)

        get_chunk_offset_shape(ndims, local_offset, shape,
                               chunk, nchunks, 
                               c, chunk_offset_in, chunk_offset_out, chunk_shape,
                               &partial_chunk);

        printf("%zu %zu\n", chunk_offset_in[0], chunk_shape[0]);
        NCERR(nc_get_vara(ncid, varid, chunk_offset_in, chunk_shape, buffer));

        MPI_Request request;
        write_uncompressed_async(var, ndims, chunk_offset_out, chunk_shape, buffer,
                                 type, async_writer_rank, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    free(buffer);
}

// Copy all chunked variables
void copy_chunked(const char * filename, int async_writer_rank) {
    int ncid;
    NCERR(nc_open(filename, NC_NOWRITE, &ncid));

    int nvars;
    NCERR(nc_inq_nvars(ncid, &nvars));

    for (int v=0; v<nvars; ++v) {
        int ndims;
        NCERR(nc_inq_varndims(ncid, v, &ndims));

        bool coll = is_collated(ncid, v);
        if (!coll) continue;

        int storage;
        size_t in_chunk[ndims];
        NCERR(nc_inq_var_chunking(ncid, v, &storage, in_chunk));
        if (storage != NC_CHUNKED) continue;

        char varname[NC_MAX_NAME+1];
        NCERR(nc_inq_varname(ncid, v, varname));

        // Get a handle to the variable on the writer
        varid_t var = open_variable_async(varname, NC_MAX_NAME+1, async_writer_rank);

        // Get the chunk info
        size_t out_chunk[ndims];
        variable_info_async(var, ndims, out_chunk, async_writer_rank);

        // Offset of this file in the collation
        size_t out_offset[ndims];
        size_t local_size[ndims];
        size_t total_size[ndims];
        {
            size_t in_offset[ndims];
            get_collation_info(ncid, v, in_offset, out_offset, local_size, total_size, ndims);
        }

        bool is_aligned = true;
        for (int d=0;d<ndims;++d) {
            is_aligned &= in_chunk[d] == out_chunk[d];
            // Start lines up with chunks
            is_aligned &= out_offset[d]%out_chunk[d] == 0;
            // End lines up with chunks, or is the final chunk
            is_aligned &= ((out_offset[d] + local_size[d])%out_chunk[d] == 0
                           || out_offset[d] + local_size[d] == total_size[d]);
        }

        if (is_aligned) {
            fprintf(stdout, "\tAligned HDF5 copy of %s from %s\n", varname, filename);
            copy_netcdf_variable_chunks(var, ncid, v, ndims, out_chunk, async_writer_rank);
        } else {
            fprintf(stdout, "\tUnaligned NetCDF4 copy of %s from %s\n", varname, filename);
            copy_netcdf_variable_chunks(var, ncid, v, ndims, out_chunk, async_writer_rank);
        }

        close_variable_async(var, async_writer_rank);
    }
}

