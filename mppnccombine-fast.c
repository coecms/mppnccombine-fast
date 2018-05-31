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
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include <mpi.h>

#include "error.h"
#include "async.h"
#include "read_chunked.h"

#define TAG_DECOMP 1
#define TAG_CHUNK 2


// Print diagnostic information about this file's collation
void print_offsets(size_t out_offset[], size_t local_size[], int ndims) {
    /*
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
    */
}

// Copy a (possibly collated) field in NetCDF mode
void copy_netcdf(int ncid_out, int varid_out, int ncid_in, int varid_in) {
    int ndims;
    NCERR(nc_inq_varndims(ncid_in, varid_in, &ndims));

    size_t in_offset[ndims];
    size_t out_offset[ndims];
    size_t local_size[ndims];
    size_t total_size[ndims];

    get_collation_info(ncid_in, varid_in, out_offset, local_size, total_size, ndims);

    size_t size = 1;
    for (int d=0; d<ndims; ++d) {
        size *= local_size[d];
        in_offset[d] = 0;
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

    int ndims;
    int nvars;
    int natts;
    NCERR(nc_inq(in_file, &ndims, &nvars, &natts, NULL));
    
    // Copy global attributes
    copy_attrs(out_file, NC_GLOBAL, in_file, NC_GLOBAL, natts);

    // Copy dimensions
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
        if (is_collated(in_file, v)) {
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


// Copy contiguous variables - no chunking means no compression, so we can just
// use NetCDF
void copy_contiguous(const char * out_path, char ** in_paths, int n_in) {

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

void check_chunking(char ** in_paths, int n_in) {
    int ncid0, nvars;
    NCERR(nc_open(in_paths[0], NC_NOWRITE, &ncid0));
    NCERR(nc_inq_nvars(ncid0, &nvars));
    
    for (int v=0; v<nvars; ++v) {
        int ndims;
        NCERR(nc_inq_varndims(ncid0, v, &ndims));

        int storage0;
        size_t chunk0[ndims];
        int shuffle0;
        int deflate0;
        int deflate_level0;
        NCERR(nc_inq_var_chunking(ncid0, v, &storage0, chunk0));
        NCERR(nc_inq_var_deflate(ncid0, v, &shuffle0, &deflate0, &deflate_level0));

        for (int i=1; i<n_in; ++i) {
            int ncid;
            NCERR(nc_open(in_paths[i], NC_NOWRITE, &ncid));
            int storage;
            size_t chunk[ndims];
            int shuffle;
            int deflate;
            int deflate_level;
            NCERR(nc_inq_var_chunking(ncid, v, &storage, chunk));
            NCERR(nc_inq_var_deflate(ncid, v, &shuffle, &deflate, &deflate_level));

            assert(storage == storage0);
            assert(shuffle == shuffle0);
            assert(deflate == deflate0);
            assert(deflate_level == deflate_level0);

            if (storage == NC_CHUNKED) {
                for (int d=0; d<ndims; ++d){
                    assert(chunk[d] == chunk0[d]);
                }
            }

            NCERR(nc_close(ncid));
        }


    }
    NCERR(nc_close(ncid0));
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

    int current_file_idx = 0;
    MPI_Win current_file_win;
    MPI_Win_create(&current_file_idx, sizeof(current_file_idx), sizeof(current_file_idx),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &current_file_win);

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
        fprintf(stderr, "ERROR: Please run with at least 2 MPI processes\n");
        exit(-1);
    }

    const char * in_path = argv[arg_index];
    const char * out_path = args.output;

    int writer_rank = 0;

    if (comm_rank == writer_rank) {
        check_chunking(argv+arg_index, argc-arg_index);

        // Copy metadata and un-collated variables
        fprintf(stdout, "\nCopying non-collated variables\n");
        init(in_path, out_path);
        // Copy contiguous variables using NetCDF
        fprintf(stdout, "\nCopying contiguous variables\n");
        copy_contiguous(out_path, argv+arg_index, argc-arg_index);

        fprintf(stdout, "\nCopying chunked variables\n");
        double t_start = MPI_Wtime();
        size_t total_size = run_async_writer(out_path);
        double t_end = MPI_Wtime();
        double total_size_gb = total_size / pow(1024,3);

        fprintf(stdout, "\nTotal compressed size %.2f GiB | %.2f GiB / sec\n",
                total_size_gb, total_size_gb/(t_end - t_start));
    } else {
        int increment = 1;
        int my_file_idx = -1;

        // Atomic post-addition of increment to current_file_idx
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, current_file_win);
        MPI_Fetch_and_op(&increment, &my_file_idx, MPI_INT, 0, 0, MPI_SUM, current_file_win);
        MPI_Win_unlock(0, current_file_win);

        while (my_file_idx < argc-arg_index) {
            // Read chunked variables using HDF5, sending data to the async_writer to be written
            copy_chunked(argv[arg_index+my_file_idx], writer_rank);

            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, current_file_win);
            MPI_Fetch_and_op(&increment, &my_file_idx, MPI_INT, 0, 0, MPI_SUM, current_file_win);
            MPI_Win_unlock(0, current_file_win);
        }
        close_async(0);
    }

    MPI_Win_free(&current_file_win);
    return MPI_Finalize();
}
