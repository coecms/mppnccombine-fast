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

//#include "netcdf.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/*
#define NCERR(x) handle_error(x, __FILE__, __LINE__)
void handle_error(int err, const char * file, int line) {
    if (err != 0) {
        const char * message = nc_strerror(err);

        fprintf(stderr, "ERROR %s:%d %s\n", file, line, message);
        exit(-1);
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
*/

void copy(const char * in_path, const char * out_path) {
    hid_t in_file = H5Fopen(in_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t out_file = H5Fopen(out_path, H5F_ACC_RDWR, H5P_DEFAULT);

    const int ndims = 4;
    int shape[4] = {240,75,540,720};
    int chunk[4] = {1,19,135,180};

    hid_t in_var = H5Dopen(in_file, "/temp", H5P_DEFAULT);
    hid_t out_var = H5Dopen(out_file, "/temp", H5P_DEFAULT);

    // Get the number of chunks, total and in each dim
    int n_chunks = 1;
    int chunk_decomp[ndims];
    for (int d=0; d<ndims; ++d) {
        chunk_decomp[d] = shape[d] / chunk[d];
        n_chunks *= chunk_decomp[d];
    }

    size_t n_buffer = 1024^3;
    void * buffer = malloc(n_buffer);
    
    for (int c=0; c<n_chunks; ++c) {
        hsize_t offset[ndims];
        int i = c;
        for (int d=3; d>=0; --d) {
            offset[d] = (i % chunk_decomp[d]) * chunk[d];
            i /= chunk_decomp[d];
        }

        hsize_t block_size;
        H5Dget_chunk_storage_size(in_var, offset, &block_size);

        if (block_size > n_buffer) {
            n_buffer = block_size;
            buffer = realloc(buffer, n_buffer);
        }

        int filter_mask = 0;
        H5DOread_chunk(in_var, H5P_DEFAULT, offset, &filter_mask, buffer);
        H5DOwrite_chunk(out_var, H5P_DEFAULT, filter_mask, offset, block_size, buffer);

        fprintf(stdout, "%d %d %d %d\n",offset[0], offset[1], offset[2], offset[3]);
    }
}

int main(int argc, char ** argv) {

    const char * in_path = "/short/v45/aek156/access-om2/archive/01deg_jra55_ryf/output243/ocean/ocean_temp_3hourly.nc.0000";
    const char * out_path = "test.nc";

    //init(in_path, out_path);
    
    copy(in_path, out_path);

    return 0;
}
