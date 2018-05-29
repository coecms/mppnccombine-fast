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
    char * temp_filename = tempnam(NULL, "mnccf");

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
                            size_t chunk_offset[],       // Offset of the chunk local data
                            size_t chunk_shape[],        // Shape of the chunk local data
                            bool * partial               // True if partial chunk
    ) {

    *partial = false;

    for (size_t i=ndims-1; i>=0; --i) {
        size_t chunkid = c % nchunks[i];

        // Start with a full sized chunk
        chunk_shape[i] = chunk[i];

        // Count from the chunk boundary to the left of the starting point
        chunk_offset[i] = floor(local_offset[i] / chunk[i])*chunk[i]
            + chunk[i] * chunkid;

        if (chunk_offset[i] < local_offset[i]) {
            // Partial chunk on the left
            *partial = true;
            chunk_shape[i] = chunk[i] - local_offset[i] - chunk_offset[i];
            chunk_offset[i] = local_offset[i];
        }

        // Total size in dimension 'i'
        size_t total_size = (chunk[i] - local_offset[i]%chunk[i])
            + chunkid * chunk[i];

        if (total_size > local_shape[i]) {
            // Partial chunk on the right
            *partial = true;
            chunk_shape[i] -= total_size - local_shape[i];
        }

        c = c / nchunks[i];
    }
}

// Read a variable `varid` from NetCDF file `ncid`, and send it to the writer
// at `writer_rank`, using uncompressed writes
void copy_netcdf_uncompressed(varid_t var, int ncid, int varid, int ndims, const size_t chunk[]) {
    size_t local_offset[ndims];
    size_t shape[ndims];
    size_t nchunks[ndims];
    size_t total_chunks = 1;
    size_t chunk_size = 1;

    int dimids[ndims];
    NCERR(nc_inq_vardimid(ncid, varid, dimids));

    for (int i=0;i<ndims;++i) {
        NCERR(nc_inq_dimlen(ncid, dimids[i], &(shape[i])));

        // Need to account for partial chunks at the start of this file's block
        nchunks[i] = ceil(((local_offset[i] % chunk[i]) + shape[i]) / chunk[i]);

        total_chunks *= nchunks[i];
        chunk_size *= chunk[i];
    }

    void * buffer = malloc(chunk_size*sizeof(double));
    nc_type type;
    NCERR(nc_inq_vartype(ncid, varid, &type));

    for (size_t c=0; c<total_chunks; ++c) {
        size_t chunk_offset[ndims]; // Offset in the data to read from
        size_t chunk_shape[ndims];  // Shape of this chunk
        bool partial_chunk;         // Is this chunk partial? (so send uncompressed)

        get_chunk_offset_shape(ndims, local_offset, shape,
                               chunk, nchunks, 
                               c, chunk_offset, chunk_shape,
                               &partial_chunk);

        NCERR(nc_get_vara(ncid, varid, chunk_offset, chunk_shape, buffer));

        write_uncompressed_async(var, ndims, chunk_offset, chunk_shape, buffer, type);
    }

    free(buffer);
}

void read_rechunk(const char * filename, const char * varname, int writer_rank) {
    int ncid;
    NCERR(nc_open(filename, NC_NOWRITE, &ncid));

    int varid;
    NCERR(nc_inq_varid(ncid, varname, &varid));

    int ndims;
    NCERR(nc_inq_varndims(ncid, varid, &ndims));

    open_variable_async(varname, strlen(varname)+1, writer_rank);

    // Stuff to get from the writer process
    size_t chunk[ndims];
    int shuffle;
    int deflate_level;
    int data_type;

}
