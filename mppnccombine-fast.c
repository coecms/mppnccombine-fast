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

#include "netcdf.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include <mpi.h>

#include "error.h"
#include "async.h"

#define TAG_DECOMP 1
#define TAG_CHUNK 2

// Get the decomposition attribute from a variable
// If this is not a decomposed variable, decompositon[] is unchanged and returns false,
// otherwise returns true
bool get_collated_dim_decomp(int ncid, const char * varname, int decomposition[4]) {
    int varid;
    int err;
    NCERR(nc_inq_varid(ncid, varname, &varid));

    err = nc_get_att(ncid, varid, "domain_decomposition", decomposition);
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
// out_offset[ndims] - The offset in the collated array of this variable
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

// Print diagnostic information about this file's collation
void print_offsets(size_t out_offset[], size_t local_size[], int ndims) {
    fprintf(stdout, "\tStart index ");
    for (int d=0; d<ndims; ++d) {
        fprintf(stdout, "% 6zu\t", out_offset[d]);
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "\tShape       ");
    for (int d=0; d<ndims; ++d) {
        fprintf(stdout, "% 6zu\t", local_size[d]);
    }
    fprintf(stdout, "\n");
}

// Copy a (possibly collated) field in NetCDF mode
void copy_netcdf(int ncid_out, int varid_out, int ncid_in, int varid_in) {
    int ndims;
    NCERR(nc_inq_varndims(ncid_in, varid_in, &ndims));

    size_t in_offset[ndims];
    size_t out_offset[ndims];
    size_t local_size[ndims];
    size_t total_size[ndims];

    get_collation_info(ncid_in, varid_in, in_offset, out_offset, local_size, total_size, ndims);

    size_t size = 1;
    for (int d=0; d<ndims; ++d) {
        size *= local_size[d];
    }

    print_offsets(out_offset, local_size, ndims);

    // Enough size for float64_t
    void * buffer = malloc(size * 8);
    NCERR(nc_get_vara(ncid_in, varid_in, in_offset, local_size, buffer));
    NCERR(nc_put_vara(ncid_out, varid_out, out_offset, local_size, buffer));
    free(buffer);
}

// Copy NetCDF attributes, either globally or for a variable
void copy_attrs(int ncid_out, int varid_out, int ncid_in, int varid_in, int natts) {
    int buffer_len = 1024*8;
    char buffer[buffer_len];

    for (int a=0; a<natts; ++a) {
        char attname[NC_MAX_NAME+1];
        nc_type atttype;
        size_t attlen;
        NCERR(nc_inq_attname(ncid_in, varid_in, a, attname));
        NCERR(nc_inq_att(ncid_in, varid_in, attname, &atttype, &attlen));
        // Check the buffer is big enough
        assert(attlen * 8 < buffer_len);
        NCERR(nc_get_att(ncid_in, varid_in, attname, buffer));
        NCERR(nc_put_att(ncid_out, varid_out, attname, atttype, attlen, buffer));
    }
}

// Copy NetCDF headers and uncollated variables from file at in_path to file at
// out_path
void init(const char * in_path, const char * out_path) {
    int in_file;
    int out_file;

    // Open both files
    NCERR(nc_open(in_path, NC_NOWRITE, &in_file));
    NCERR(nc_create(out_path, NC_NETCDF4 | NC_CLOBBER, &out_file));

    // Copy dimensions
    int ndims;
    NCERR(nc_inq_ndims(in_file, &ndims));
    for (int d=0; d<ndims;++d) {
        char name[NC_MAX_NAME+1];
        size_t len;
        int dimid;
        int varid;

        NCERR(nc_inq_dim(in_file, d, name, &len));

        // Check if the variable with the same name is collated
        NCERR(nc_inq_varid(in_file, name, &varid));
        if (is_collated(in_file, varid)) {
            // If so get the full length
            get_collated_dim_len(in_file, name, &len);
        }

        // Create the out dim
        NCERR(nc_def_dim(out_file, name, len, &dimid));
        assert(dimid == d);
    }

    // Copy variables
    int nvars;
    NCERR(nc_inq_nvars(in_file, &nvars));
    for (int v=0; v<nvars; ++v) {
        char name[NC_MAX_NAME+1];
        nc_type type;
        int ndims;
        int natts;

        NCERR(nc_inq_var(in_file, v, name, &type, &ndims, NULL, &natts));

        int dimids[ndims];
        NCERR(nc_inq_vardimid(in_file, v, dimids));

        int out_v;
        NCERR(nc_def_var(out_file, name, type, ndims, dimids, &out_v));

        // Chunking needs to be identical in 'in' and 'out' files
        int storage;
        size_t chunk[ndims];
        NCERR(nc_inq_var_chunking(in_file, v, &storage, chunk));
        if (storage == NC_CHUNKED) {
            NCERR(nc_def_var_chunking(out_file, out_v, storage, chunk));
        }

        // Compression needs to be identical in 'in' and 'out' files
        int shuffle;
        int deflate;
        int deflate_level;
        NCERR(nc_inq_var_deflate(in_file, v, &shuffle, &deflate, &deflate_level));
        if (shuffle || deflate) {
            NCERR(nc_def_var_deflate(out_file, out_v, shuffle, deflate, deflate_level));
        }

        // Copy attributes
        copy_attrs(out_file, v, in_file, v, natts);

        // If the field is not collated copy it now
        if (! is_collated(in_file, v)) {
            fprintf(stdout, "\tUncollated NetCDF copy of %s\n", name);
            copy_netcdf(out_file, out_v, in_file, v);
        }
    }

    // Close the files
    NCERR(nc_close(in_file));
    NCERR(nc_close(out_file));
}


// Copy a variable using HDF5's optimised IO
// This is faster than the normal IO as it doesn't need to de-compress and
// re-compress the data, however it only works when the source and target
// chunking and compression settings are identical
size_t hdf5_raw_copy(
                  const size_t out_offset[], // Output offset [ndims]
                  hid_t in_var,              // Input hdf5 variable
                  const size_t in_offset[],  // Input offset [ndims]
                  const size_t shape[],      // Shape to copy [ndims]
                  int ndims                  // Number of dimensions
                 ) {
    size_t total_copied_size = 0;

    // Get the chunk metadata
    hsize_t chunk[ndims];
    hid_t in_plist = H5Dget_create_plist(in_var);
    H5ERR(H5Pget_chunk(in_plist, ndims, chunk));

    // Get the number of chunks, total and in each dim
    int n_chunks = 1;
    int chunk_decomp[ndims];
    for (int d=0; d<ndims; ++d) {
        chunk_decomp[d] = shape[d] / chunk[d];
        n_chunks *= chunk_decomp[d];
    }

    // Buffer 
    size_t n_buffer = 1024^3;
    void * buffer = malloc(n_buffer);

    hsize_t copy_out_offset[ndims];
    hsize_t copy_in_offset[ndims];

    // Loop over all the chunks
    for (int c=0; c<n_chunks; ++c) {
        hsize_t offset[ndims];
        int i = c;
        for (int d=ndims-1; d>=0; --d) {
            offset[d] = (i % chunk_decomp[d]) * chunk[d];
            i /= chunk_decomp[d];

            copy_out_offset[d] = out_offset[d] + offset[d];
            copy_in_offset[d]  = in_offset[d]  + offset[d];
        }

        // Get the block size
        hsize_t block_size;
        H5ERR(H5Dget_chunk_storage_size(in_var, copy_in_offset, &block_size));

        // Make sure the buffer is large enough
        if (block_size > n_buffer) {
            n_buffer = block_size;
            buffer = realloc(buffer, n_buffer);
        }

        // Copy this chunk's data
        uint32_t filter_mask = 0;
        H5ERR(H5DOread_chunk(in_var, H5P_DEFAULT, copy_in_offset, &filter_mask, buffer));

        MPI_Request request;
        write_chunk_async(ndims, filter_mask, copy_out_offset, block_size, buffer, 0, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);

        total_copied_size += block_size;
    }

    free(buffer);

    return total_copied_size;
}

// Copy chunked variables from the file at in_path to HDF5 variable out_var
size_t copy_chunked_variable(const char * in_path, const char * varname) {
    // Open in NetCDF mode to gather metadata
    int in_nc4;
    NCERR(nc_open(in_path, NC_NOWRITE, &in_nc4));

    int varid;
    int ndims;
    NCERR(nc_inq_varid(in_nc4, varname, &varid));
    NCERR(nc_inq_varndims(in_nc4, varid, &ndims));

    size_t in_offset[ndims];
    size_t out_offset[ndims];
    size_t local_shape[ndims];
    size_t global_shape[ndims];

    get_collation_info(in_nc4, varid, in_offset, out_offset, local_shape, global_shape, ndims);
    NCERR(nc_close(in_nc4));

    fprintf(stdout, "\tHDF5 copy of %s from %s\n", varname, in_path);
    print_offsets(out_offset, local_shape, ndims);

    // Open in HDF5 mode to do the copy
    hid_t in_file = H5Fopen(in_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    H5ERR(in_file);
    hid_t in_var = H5Dopen(in_file, varname, H5P_DEFAULT);
    H5ERR(in_var);

    size_t total_copied_size = hdf5_raw_copy(out_offset, in_var, in_offset, local_shape, ndims);

    H5ERR(H5Dclose(in_var));
    H5ERR(H5Fclose(in_file));

    return total_copied_size;
}

// Copy chunked variables - these may be compressed, so we'll use HDF5. Since
// we can't have the same file open in both HDF5 and NetCDF4 modes we need to
// do a bit of shuffling to get all the metadata.
void copy_chunked(char ** in_paths, int n_in) {
    size_t total_copied_size = 0;
    double t_start = MPI_Wtime();

    // Get the total number of variables
    int in_nc4;
    NCERR(nc_open(in_paths[0], NC_NOWRITE, &in_nc4));
    int nvars;
    NCERR(nc_inq_nvars(in_nc4, &nvars));
    NCERR(nc_close(in_nc4));

    fprintf(stdout, "\nCopying chunked variables\n");

    // Loop over each variable
    for (int v=0; v<nvars; ++v) {
        int in_nc4;
        NCERR(nc_open(in_paths[0], NC_NOWRITE, &in_nc4));
        char varname[NC_MAX_NAME+1];
        NCERR(nc_inq_varname(in_nc4, v, varname));
        int storage;
        NCERR(nc_inq_var_chunking(in_nc4, v, &storage, NULL));
        bool coll = is_collated(in_nc4, v);
        NCERR(nc_close(in_nc4));

        if (!coll) continue;

        // Copy chunked variable from all input files
        if (storage == NC_CHUNKED) {
            open_variable_async(varname, NC_MAX_NAME+1, 0);
            for (int i=0; i<n_in; ++i) {
                total_copied_size += copy_chunked_variable(in_paths[i], varname);
            }
            close_variable_async(0);
        }
    }

    double t_end = MPI_Wtime();
    double total_size_gb = total_copied_size / pow(1024,3);

    fprintf(stdout, "Total compressed size %.2f GiB | %.2f GiB / sec\n", total_size_gb, total_size_gb/(t_end - t_start));
}

// Copy contiguous variables - no chunking means no compression, so we can just
// use NetCDF
void copy_contiguous(const char * out_path, char ** in_paths, int n_in) {
    fprintf(stdout, "\nCopying contiguous variables\n");

    int out_nc4;
    NCERR(nc_open(out_path, NC_WRITE, &out_nc4));

    int nvars;
    NCERR(nc_inq_nvars(out_nc4, &nvars));

    for (int v=0; v<nvars; ++v) {
        // Check the metadata of the output file to see if this variable is contiguous
        char varname[NC_MAX_NAME+1];
        NCERR(nc_inq_varname(out_nc4, v, varname));
        int storage;
        NCERR(nc_inq_var_chunking(out_nc4, v, &storage, NULL));
        if (storage == NC_CONTIGUOUS) {
            if (is_collated(out_nc4, v)) {
                for (int i=0; i<n_in; ++i) {
                    fprintf(stdout, "\tNetCDF copy of %s from %s\n", varname, in_paths[i]);
                    int in_nc4;
                    NCERR(nc_open(in_paths[i], NC_NOWRITE, &in_nc4));
                    copy_netcdf(out_nc4, v, in_nc4, v);
                    NCERR(nc_close(in_nc4));
                }
            }
        }
    }
    NCERR(nc_close(out_nc4));
}

// Consume the variable from a single input file
size_t copy_mpi_consume_file_variable(hid_t out_var) {
    size_t total_copied_size = 0;
    int ndims;
    MPI_Recv(&ndims, 1, MPI_INT, 1, TAG_DECOMP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    uint64_t out_offset[ndims];
    uint64_t local_shape[ndims];
    MPI_Recv(out_offset, ndims, MPI_UINT64_T, 1, TAG_DECOMP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(local_shape, ndims, MPI_UINT64_T, 1, TAG_DECOMP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Get the chunk metadata
    hsize_t chunk[ndims];
    hid_t out_plist = H5Dget_create_plist(out_var);
    H5ERR(H5Pget_chunk(out_plist, ndims, chunk));

    // Get the number of chunks, total and in each dim
    size_t n_chunks = 1;
    int chunk_decomp[ndims];
    for (int d=0; d<ndims; ++d) {
        chunk_decomp[d] = local_shape[d] / chunk[d];
        n_chunks *= chunk_decomp[d];
    }

    // Buffer 
    uint64_t n_buffer = 1024^3;
    void * buffer = malloc(n_buffer);

    hsize_t copy_out_offset[ndims];

    // Loop over all the chunks
    for (size_t c=0; c<n_chunks; ++c) {
        hsize_t offset[ndims];
        int i = c;
        for (int d=ndims-1; d>=0; --d) {
            offset[d] = (i % chunk_decomp[d]) * chunk[d];
            i /= chunk_decomp[d];

            copy_out_offset[d]  = out_offset[d]  + offset[d];
        }

        // Get the block size
        uint64_t block_size;
        MPI_Recv(&block_size, 1, MPI_UINT64_T, 1, TAG_CHUNK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        uint32_t filter_mask = 0;
        MPI_Recv(&filter_mask, 1, MPI_UINT32_T, 1, TAG_CHUNK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


        // Make sure the buffer is large enough
        if (block_size > n_buffer) {
            n_buffer = block_size;
            buffer = realloc(buffer, n_buffer);
        }

        MPI_Recv(buffer, block_size, MPI_CHAR, 1, TAG_CHUNK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Copy this chunk's data
        H5ERR(H5DOwrite_chunk(out_var, H5P_DEFAULT, filter_mask, copy_out_offset, block_size, buffer));

        total_copied_size += block_size;
    }

    free(buffer);

    MPI_Barrier(MPI_COMM_WORLD);

    return total_copied_size;
}

// Produce the variable from a single input file
void copy_mpi_produce_file_variable(const char * in_path, const  char * varname) {
    fprintf(stdout, "\tHDF5 copy of %s from %s\n", varname, in_path);
    // Open in NetCDF mode to gather metadata
    int in_nc4;
    NCERR(nc_open(in_path, NC_NOWRITE, &in_nc4));

    int varid;
    int ndims;
    NCERR(nc_inq_varid(in_nc4, varname, &varid));
    NCERR(nc_inq_varndims(in_nc4, varid, &ndims));

    size_t in_offset[ndims];
    uint64_t out_offset[ndims];
    uint64_t local_shape[ndims];
    uint64_t global_shape[ndims];

    get_collation_info(in_nc4, varid, in_offset, out_offset, local_shape, global_shape, ndims);
    NCERR(nc_close(in_nc4));

    // Send the out_offset to the consumer
    MPI_Send(&ndims, 1, MPI_INT, 0, TAG_DECOMP, MPI_COMM_WORLD);
    MPI_Send(out_offset, ndims, MPI_UINT64_T, 0, TAG_DECOMP, MPI_COMM_WORLD);
    MPI_Send(local_shape, ndims, MPI_UINT64_T, 0, TAG_DECOMP, MPI_COMM_WORLD);

    // Open in HDF5 mode to do the copy
    hid_t in_file = H5Fopen(in_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    H5ERR(in_file);
    hid_t in_var = H5Dopen(in_file, varname, H5P_DEFAULT);
    H5ERR(in_var);

    // Get the chunk metadata
    hsize_t chunk[ndims];
    hid_t in_plist = H5Dget_create_plist(in_var);
    H5ERR(H5Pget_chunk(in_plist, ndims, chunk));

    // Get the number of chunks, total and in each dim
    int n_chunks = 1;
    int chunk_decomp[ndims];
    for (int d=0; d<ndims; ++d) {
        chunk_decomp[d] = local_shape[d] / chunk[d];
        n_chunks *= chunk_decomp[d];
    }

    // Double buffer so we can use Isend
    uint64_t n_buffer[2] = {0};
    void * buffer[2] = {0};
    uint32_t filter_mask[2];
    uint64_t block_size_64[2];
    MPI_Request request[2][3];

    hsize_t copy_in_offset[ndims];

    // Loop over all the chunks
    for (int c=0; c<n_chunks; ++c) {
        hsize_t offset[ndims];
        int i = c;
        for (int d=ndims-1; d>=0; --d) {
            offset[d] = (i % chunk_decomp[d]) * chunk[d];
            i /= chunk_decomp[d];

            copy_in_offset[d]  = in_offset[d]  + offset[d];
        }

        // Get the block size
        hsize_t block_size;
        H5ERR(H5Dget_chunk_storage_size(in_var, copy_in_offset, &block_size));

        if (c >= 2) {
            // Wait for buffer to catch up
            MPI_Waitall(3, request[c%2], MPI_STATUS_IGNORE);
        }

        // Make sure the buffer is large enough
        if (block_size > n_buffer[c%2]) {
            n_buffer[c%2] = block_size;
            buffer[c%2] = realloc(buffer[c%2], n_buffer[c%2]);
        }

        // Copy this chunk's data
        H5ERR(H5DOread_chunk(in_var, H5P_DEFAULT, copy_in_offset, &filter_mask[c%2], buffer[c%2]));
        
        // Send
        block_size_64[c%2] = block_size;
        MPI_Isend(&block_size_64[c%2], 1, MPI_UINT64_T, 0, TAG_CHUNK, MPI_COMM_WORLD, &(request[c%2][0]));
        MPI_Isend(&filter_mask[c%2], 1, MPI_UINT32_T, 0, TAG_CHUNK, MPI_COMM_WORLD, &(request[c%2][1]));
        MPI_Isend(buffer[c%2], block_size, MPI_CHAR, 0, TAG_CHUNK, MPI_COMM_WORLD, &(request[c%2][2]));
    }

    free(buffer[0]);
    free(buffer[1]);
    H5ERR(H5Dclose(in_var));
    H5ERR(H5Fclose(in_file));
    MPI_Barrier(MPI_COMM_WORLD);
}

void copy_chunked_mpi_consumer(const char * out_path, size_t n_in) {
    size_t total_copied_size = 0;
    double t_start = MPI_Wtime();

    // Get the total number of variables
    int out_nc4;
    NCERR(nc_open(out_path, NC_NOWRITE, &out_nc4));
    int nvars;
    NCERR(nc_inq_nvars(out_nc4, &nvars));
    NCERR(nc_close(out_nc4));

    fprintf(stdout, "\nCopying chunked variables\n");

    char varname[NC_MAX_NAME+1];

    // Loop over each variable
    for (int v=0; v<nvars; ++v) {
        int out_nc4;
        NCERR(nc_open(out_path, NC_NOWRITE, &out_nc4));
        NCERR(nc_inq_varname(out_nc4, v, varname));
        int storage;
        NCERR(nc_inq_var_chunking(out_nc4, v, &storage, NULL));
        bool coll = is_collated(out_nc4, v);
        NCERR(nc_close(out_nc4));

        if (!coll) continue;

        // Copy chunked variable from all input files
        if (storage == NC_CHUNKED) {
            // Open the output file in HDF5 mode
            hid_t out_file = H5Fopen(out_path, H5F_ACC_RDWR, H5P_DEFAULT);
            H5ERR(out_file);

            // Request the variable from the producer
            MPI_Bcast(varname, NC_MAX_NAME+1, MPI_CHAR, 0, MPI_COMM_WORLD);

            hid_t out_var = H5Dopen(out_file, varname, H5P_DEFAULT);
            H5ERR(out_var);
            for (int i=0; i<n_in; ++i) {
                total_copied_size += copy_mpi_consume_file_variable(out_var);
            }
            H5ERR(H5Dclose(out_var));
            H5ERR(H5Fclose(out_file));
        }
    }

    // Tell the producer we're finished
    varname[0] = '\0';
    MPI_Bcast(varname, NC_MAX_NAME+1, MPI_CHAR, 0, MPI_COMM_WORLD);

    double t_end = MPI_Wtime();
    double total_size_gb = total_copied_size / pow(1024,3);

    fprintf(stdout, "Total compressed size %.2f GiB | %.2f GiB / sec\n", total_size_gb, total_size_gb/(t_end - t_start));
}

void copy_chunked_mpi_producer(char ** in_paths, size_t n_in) {
    while (true) {
        // Variable requested by consumer
        char varname[NC_MAX_NAME+1];
        MPI_Bcast(varname, NC_MAX_NAME+1, MPI_CHAR, 0, MPI_COMM_WORLD);

        if (strlen(varname) == 0) {
            // No more requests
            return;
        }
        
        for (int i=0; i<n_in; ++i) {
            // Send each file's variable to the consumer
            copy_mpi_produce_file_variable(in_paths[i], varname);
        }
    }
}

// Split into a producer that reads the input files and a consumer that writes
// to the output file
void copy_chunked_mpi(const char * out_path, char ** in_paths, size_t n_in, int comm_rank) {
    if (comm_rank == 0) {
        copy_chunked_mpi_consumer(out_path, n_in);
    } else {
        copy_chunked_mpi_producer(in_paths, n_in);
    }
}

struct args_t {
    const char * output;
};

static char doc[] = "Quickly collate MOM output files";

static struct argp_option opts[] = {
    {"output", 'o', "FILE", 0, "Output file"},
    {0},
};

static error_t parse_opt(int key, char *arg, struct argp_state * state) {
    struct args_t * args = state->input;

    switch(key) {
        case 'o':
            args->output = arg;
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static struct argp argp = {
    .parser = parse_opt,
    .options = opts,
    .doc = doc,
};

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    int comm_rank;
    int comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int arg_index;
    struct args_t args = {0};

    argp_parse(&argp, argc, argv, 0, &arg_index, &args);
    if (args.output == NULL) {
        fprintf(stderr, "ERROR: No output file specified\n");
        exit(-1);
    }
    if (arg_index == argc) {
        fprintf(stderr, "ERROR: No input files specified\n");
        exit(-1);
    }
    if (comm_size < 2) {
        fprintf(stderr, "ERROR: Please run with `mpirun -n 2`\n");
        exit(-1);
    }

    const char * in_path = argv[arg_index];
    const char * out_path = args.output;

    if (comm_rank == 0) {
        // Copy metadata and un-collated variables
        init(in_path, out_path);
        // Copy contiguous variables using NetCDF
        copy_contiguous(out_path, argv+arg_index, argc-arg_index);

        run_async_writer(out_path);
    } else {

        // Copy chunked variables using HDF5
        copy_chunked(argv+arg_index+comm_rank-1, 1);
        close_async(0);
    }

    return MPI_Finalize();
}
