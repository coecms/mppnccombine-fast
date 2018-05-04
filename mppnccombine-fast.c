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

#define NCERR(x) handle_nc_error(x, __FILE__, __LINE__)
void handle_nc_error(int err, const char * file, int line) {
    if (err != 0) {
        const char * message = nc_strerror(err);

        fprintf(stderr, "ERROR %s:%d %d %s\n", file, line, err, message);
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

        int decomposition[4];
        is_collated[d] = get_collated_dim_decomp(ncid, dimname, decomposition);

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

        out = out || is_collated[d];
    }

    return out;
}

// Returns true if any of the dimensions are collated
bool is_collated(int ncid, int varid) {
    int ndims;
    NCERR(nc_inq_varndims(ncid, varid, &ndims));

    size_t in_offset[ndims];
    size_t out_offset[ndims];
    size_t local_size[ndims];
    size_t total_size[ndims];

    return get_collation_info(ncid, varid, in_offset, out_offset, local_size, total_size, ndims);
}

// Copy an uncollated field
/*
void copy_netcdf(int ncid_out, int varid_out, int ncid_in, int varid_in) {
    int ndims;
    NCERR(nc_inq_varndims(ncid_in, varid_in, &ndims));

    int dimids[ndims];
    NCERR(nc_inq_vardimid(ncid_in, varid_in, dimids));

    size_t size = 1;
    for (int d=0; d<ndims; ++d) {
        size_t len;
        NCERR(nc_inq_dimlen(ncid_in, dimids[d], &len));
        size *= len;
    }

    void * buffer = malloc(size * 8);
    NCERR(nc_get_var(ncid_in, varid_in, buffer));
    NCERR(nc_put_var(ncid_out, varid_out, buffer));
    free(buffer);
}
*/

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

    void * buffer = malloc(size * 8);
    NCERR(nc_get_vara(ncid_in, varid_in, in_offset, local_size, buffer));
    NCERR(nc_put_vara(ncid_out, varid_out, out_offset, local_size, buffer));
    free(buffer);
}

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

        int varid;
        NCERR(nc_inq_varid(in_file, name, &varid));

        if (is_collated(in_file, varid)) {
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
        if (storage == NC_CHUNKED) {
            fprintf(stdout, "%s chunked\n", name);
            NCERR(nc_def_var_chunking(out_file, out_v, storage, chunk));
        }

        int shuffle;
        int deflate;
        int deflate_level;
        NCERR(nc_inq_var_deflate(in_file, v, &shuffle, &deflate, &deflate_level));
        if (shuffle || deflate) {
            NCERR(nc_def_var_deflate(out_file, out_v, shuffle, deflate, deflate_level));
        }

        int no_fill;
        double fill_buffer[1] = {0,};
        NCERR(nc_inq_var_fill(in_file, v, &no_fill, fill_buffer));
        if (!no_fill) {
            NCERR(nc_def_var_fill(out_file, out_v, no_fill, fill_buffer));
        }
        fprintf(stdout, "%s %d %f\n", name, no_fill, fill_buffer[0]);

        copy_attrs(out_file, v, in_file, v, natts);

        if (! is_collated(in_file, v)) {
            fprintf(stdout, "NetCDF copy of %s\n", name);
            copy_netcdf(out_file, out_v, in_file, v);
        }
    }


    NCERR(nc_close(in_file));
    NCERR(nc_close(out_file));
}


int hdf5_raw_copy(
                  hid_t out_var,          // Output hdf5 variable
                  const size_t out_offset[], // Output offset [ndims]
                  hid_t in_var,           // Input hdf5 variable
                  const size_t in_offset[],  // Input offset [ndims]
                  const size_t shape[],      // Shape to copy [ndims]
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
        for (int d=ndims-1; d>=0; --d) {
            offset[d] = (i % chunk_decomp[d]) * chunk[d];
            i /= chunk_decomp[d];

            copy_out_offset[d] = out_offset[d] + offset[d];
            copy_in_offset[d]  = in_offset[d]  + offset[d];
        }

        /*
        // Debugging
        if (c < 2) {
        fprintf(stdout, "in  % 4zu % 4zu % 4zu % 4zu\n", copy_in_offset[0], copy_in_offset[1], copy_in_offset[2], copy_in_offset[3]);
        fprintf(stdout, "out % 4zu % 4zu % 4zu % 4zu\n", copy_out_offset[0], copy_out_offset[1], copy_out_offset[2], copy_out_offset[3]);
        }
        */

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
    }

    return 0;
}

void chunked_copy(const char * in_path, hid_t out_var, const char * varname) {
    int in_nc4;
    NCERR(nc_open(in_path, NC_NOWRITE, &in_nc4));

    int varid;
    NCERR(nc_inq_varid(in_nc4, varname, &varid));

    int ndims;
    NCERR(nc_inq_varndims(in_nc4, varid, &ndims));

    size_t in_offset[ndims];
    size_t out_offset[ndims];
    size_t local_shape[ndims];
    size_t global_shape[ndims];

    get_collation_info(in_nc4, varid, in_offset, out_offset, local_shape, global_shape, ndims);
    NCERR(nc_close(in_nc4));

    fprintf(stdout, "HDF5 copy of %s from %s\n", varname, in_path);
    fprintf(stdout, "\tStart index ");
    for (int d=0; d<ndims; ++d) {
        fprintf(stdout, "% 6zu\t", out_offset[d]);
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "\tShape       ");
    for (int d=0; d<ndims; ++d) {
        fprintf(stdout, "% 6zu\t", local_shape[d]);
    }
    fprintf(stdout, "\n");


    hid_t in_file = H5Fopen(in_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    H5ERR(in_file);

    hid_t in_var = H5Dopen(in_file, varname, H5P_DEFAULT);
    H5ERR(in_var);

    hdf5_raw_copy(out_var, out_offset, in_var, in_offset, local_shape, ndims);

    H5ERR(H5Fflush(out_var, H5F_SCOPE_GLOBAL));

    H5ERR(H5Dclose(in_var));
    H5ERR(H5Fclose(in_file));
}

void copy_contiguous(const char * out_path, char ** in_paths, int n_in) {
    fprintf(stdout, "Copying contiguous variables\n");

    int out_nc4;
    NCERR(nc_open(out_path, NC_WRITE, &out_nc4));

    int nvars;
    NCERR(nc_inq_nvars(out_nc4, &nvars));

    for (int v=0; v<nvars; ++v) {
        char varname[NC_MAX_NAME+1];
        NCERR(nc_inq_varname(out_nc4, v, varname));
        int storage;
        NCERR(nc_inq_var_chunking(out_nc4, v, &storage, NULL));
        if (storage == NC_CONTIGUOUS) {
            if (is_collated(out_nc4, v)) {
                for (int i=0; i<n_in; ++i) {
                    fprintf(stdout, "Copying %s from %s\n", varname, in_paths[i]);
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

    init(in_path, out_path);

    int in_nc4;
    NCERR(nc_open(in_path, NC_NOWRITE, &in_nc4));
    int nvars;
    NCERR(nc_inq_nvars(in_nc4, &nvars));
    NCERR(nc_close(in_nc4));

    copy_contiguous(out_path, argv+arg_index, argc-arg_index);

    fprintf(stdout, "Copying chunked variables\n");

    hid_t out_file = H5Fopen(out_path, H5F_ACC_RDWR, H5P_DEFAULT);
    H5ERR(out_file);

    for (int v=0; v<nvars; ++v) {
        int in_nc4;
        NCERR(nc_open(in_path, NC_NOWRITE, &in_nc4));
        char varname[NC_MAX_NAME+1];
        NCERR(nc_inq_varname(in_nc4, v, varname));
        int storage;
        NCERR(nc_inq_var_chunking(in_nc4, v, &storage, NULL));
        bool coll = is_collated(in_nc4, v);
        NCERR(nc_close(in_nc4));

        if (!coll) continue;

        // Copy chunked variables
        if (storage == NC_CHUNKED) {
            hid_t out_var = H5Dopen(out_file, varname, H5P_DEFAULT);
            H5ERR(out_var);
            for (int i=arg_index; i<argc; ++i) {
                chunked_copy(argv[i], out_var, varname);
            }
            H5ERR(H5Dclose(out_var));
        }
    }
    H5ERR(H5Fclose(out_file));

    return 0;
}
