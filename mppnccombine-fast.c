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
#include <mpi.h>

#define NCERR(x) handle_nc_error(x, __FILE__, __LINE__)
void handle_nc_error(int err, const char * file, int line) {
    if (err != 0) {
        const char * message = nc_strerror(err);

        fprintf(stderr, "ERROR %s:%d %s\n", file, line, message);
        exit(-1);
    }
}


#define H5ERR(x) handle_h5_error(x, __FILE__, __LINE__)
void handle_h5_error(int err, const char * file, int line) {
    if (err < 0) {
        fprintf(stderr, "ERROR %s:%d\n", file, line);
        H5Eprint1(stderr);
        exit(-1);
    }
}

void get_collated_dim_decomp(int ncid, const char * varname, int decomposition[4]) {
    int varid;
    NCERR(nc_inq_varid(ncid, varname, &varid)); 
    NCERR(nc_get_att(ncid, varid, "domain_decomposition", decomposition));
}

void get_collated_dim_len(int ncid, const char * varname, size_t * len) {
    int decomposition[4];
    get_collated_dim_decomp(ncid, varname, decomposition);
    *len = decomposition[1];
}

// Get the output offset for a 4d variable
void get_out_offset_4d(int ncid, int out_offset[4]) {
    out_offset[0] = 0;
    out_offset[1] = 0;
    
    int decomposition[4];
    get_collated_dim_decomp(ncid, "yt_ocean", decomposition);
    out_offset[2] = decomposition[2]-1;

    get_collated_dim_decomp(ncid, "xt_ocean", decomposition);
    out_offset[3] = decomposition[2]-1;
}

void init(const char * in_path, const char * out_path) {
    int in_file;
    int out_file;

    NCERR(nc_open(in_path, NC_NOWRITE, &in_file));
    NCERR(nc_create(out_path, NC_NETCDF4 | NC_CLOBBER, &out_file));

    // Copy dimensions
    int ndims;
    NCERR(nc_inq_ndims(in_file, &ndims));
    for (int d=0; d<ndims;++d) {
        char name[NC_MAX_NAME+1];
        size_t len;
        int dimid;

        NCERR(nc_inq_dim(in_file, d, name, &len));

        if (strcmp(name, "xt_ocean") == 0) {
            get_collated_dim_len(in_file, name, &len);
        }
        if (strcmp(name, "yt_ocean") == 0) {
            get_collated_dim_len(in_file, name, &len);
        }

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

        int storage;
        size_t chunk[ndims];
        NCERR(nc_inq_var_chunking(in_file, v, &storage, chunk));
        NCERR(nc_def_var_chunking(out_file, out_v, storage, chunk));

        int shuffle;
        int deflate;
        int deflate_level;
        NCERR(nc_inq_var_deflate(in_file, v, &shuffle, &deflate, &deflate_level));
        NCERR(nc_def_var_deflate(out_file, out_v, shuffle, deflate, deflate_level));

        int no_fill;
        double fill_buffer[1] = {0,};
        NCERR(nc_inq_var_fill(in_file, v, &no_fill, fill_buffer));
        NCERR(nc_def_var_fill(out_file, out_v, no_fill, fill_buffer));
        fprintf(stdout, "%s %f\n", name, fill_buffer[0]);
    }


    NCERR(nc_close(in_file));
    NCERR(nc_close(out_file));
}


int hdf5_raw_copy(
        hid_t out_var,          // Output hdf5 variable
        const int out_offset[], // Output offset [ndims]
        hid_t in_var,           // Input hdf5 variable
        const int in_offset[],  // Input offset [ndims]
        const int shape[],      // Shape to copy [ndims]
        int ndims               // Number of dimensions
) {
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
        for (int d=3; d>=0; --d) {
            offset[d] = (i % chunk_decomp[d]) * chunk[d];
            i /= chunk_decomp[d];

            copy_out_offset[d] = out_offset[d] + offset[d];
            copy_in_offset[d]  = in_offset[d]  + offset[d];
        }

        if (c < 2) {
            fprintf(stdout, "in  % 4zu % 4zu % 4zu % 4zu\n", copy_in_offset[0], copy_in_offset[1], copy_in_offset[2], copy_in_offset[3]);
            fprintf(stdout, "out % 4zu % 4zu % 4zu % 4zu\n", copy_out_offset[0], copy_out_offset[1], copy_out_offset[2], copy_out_offset[3]);
        }

        // Get the block size
        hsize_t block_size;
        H5ERR(H5Dget_chunk_storage_size(in_var, copy_in_offset, &block_size));

        // Make sure the buffer is large enough
        if (block_size > n_buffer) {
            n_buffer = block_size;
            buffer = realloc(buffer, n_buffer);
        }

        // Copy this chunk's block
        uint32_t filter_mask = 0;
        H5ERR(H5DOread_chunk(in_var, H5P_DEFAULT, copy_in_offset, &filter_mask, buffer));
        H5ERR(H5DOwrite_chunk(out_var, H5P_DEFAULT, filter_mask, copy_out_offset, block_size, buffer));

        if (c % 1000 == 99) {
            fprintf(stdout, "Flushing at %d/%d\n", c, n_chunks);
            //H5Dflush(out_var);
            break;
        }
    }

    return 0;
}

void copy(const char * in_path, hid_t out_var) {
    int err;
    int in_nc4;
    int out_offset[4] = {0,0,0,0};


    NCERR(nc_open(in_path, NC_NOWRITE, &in_nc4));
    get_out_offset_4d(in_nc4, out_offset);
    NCERR(nc_close(in_nc4));

    fprintf(stdout, "% 4d % 4d % 4d % 4d %s\n", out_offset[0], out_offset[1], out_offset[2], out_offset[3], in_path);

    hid_t in_file = H5Fopen(in_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    H5ERR(in_file);

    const int ndims = 4;
    int shape[4] = {240,75,540,720};

    hid_t in_var = H5Dopen(in_file, "/temp", H5P_DEFAULT);
    H5ERR(in_var);

    int in_offset[4] = {0,0,0,0};

    hdf5_raw_copy(out_var, out_offset, in_var, in_offset, shape, ndims);

    H5ERR(H5Dclose(in_var));

    H5ERR(H5Fclose(in_file));
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

    int comm_size;
    int comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

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

    const char * in_path = argv[arg_index];
    const char * out_path = args.output;

    if (comm_rank == 0) {
        // Setup the headers and non-shared data fields
        init(in_path, out_path);
    }

    hid_t out_file_props = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(out_file_props, MPI_COMM_WORLD, MPI_INFO_NULL);

    hid_t out_file = H5Fopen(out_path, H5F_ACC_RDWR, out_file_props);
    H5ERR(out_file);

    hid_t out_var = H5Dopen(out_file, "/temp", H5P_DEFAULT);
    H5ERR(out_var);

    fprintf(stdout, "Rank %d has opened\n", comm_rank);
    
    // Loop over each input file, copying shared data fields
    for (int i=arg_index + comm_rank; i<argc; i+=comm_size) {
        copy(argv[i], out_var);
    }

    fprintf(stdout, "Rank %d is closing\n", comm_rank);

    H5ERR(H5Dclose(out_var));
    H5ERR(H5Fclose(out_file));

    MPI_Finalize();

    return 0;
}
