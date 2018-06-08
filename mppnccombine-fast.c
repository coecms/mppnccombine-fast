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


// Returns true if dimension dimid in file ncid is unlimited
bool is_unlimited(int ncid, int dimid) {
    int nudims;
    NCERR(nc_inq_unlimdims(ncid, &nudims, NULL));
    int udims[nudims];
    NCERR(nc_inq_unlimdims(ncid, NULL, udims));

    for (int i=0; i<nudims; ++i) {
        if (udims[i] == dimid) {
            return true;
        }
    }
    return false;
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
        if (varid_in == NC_GLOBAL && strcmp(attname, "NumFilesInSet") == 0)
            continue; // Don't copy
        // Check the buffer is big enough
        assert(attlen * 8 < buffer_len);
        NCERR(nc_get_att(ncid_in, varid_in, attname, buffer));
        if (varid_in == NC_GLOBAL && strcmp(attname, "filename") == 0)
        {
            // Replace filename attribute with new collated filename
            NCERR(nc_inq_path(ncid_out, &attlen, NULL));
            // Check the buffer is big enough
            assert(attlen * 8 < buffer_len);
            NCERR(nc_inq_path(ncid_out, NULL, buffer));
        }
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


// Copy a variable using HDF5's optimised IO
// This is faster than the normal IO as it doesn't need to de-compress and
// re-compress the data, however it only works when the source and target
// chunking and compression settings are identical
size_t hdf5_raw_copy(
                  varid_t varid,
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
        write_chunk_async(varid, ndims, filter_mask, copy_out_offset, block_size, buffer, 0, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);

        total_copied_size += block_size;
    }

    free(buffer);

    return total_copied_size;
}

// Copy chunked variables from the file at in_path to HDF5 variable out_var
size_t copy_chunked_variable(varid_t var_id, const char * in_path, const char * varname) {
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

    size_t total_copied_size = hdf5_raw_copy(var_id, out_offset, in_var, in_offset, local_shape, ndims);

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
            varid_t var_id = open_variable_async(varname, NC_MAX_NAME+1, 0);
            for (int i=0; i<n_in; ++i) {
                total_copied_size += copy_chunked_variable(var_id, in_paths[i], varname);
            }
            close_variable_async(var_id, 0);
        }
    }

    double t_end = MPI_Wtime();
    double total_size_gb = total_copied_size / pow(1024,3);

    // fprintf(stdout, "Total compressed size %.2f GiB | %.2f GiB / sec\n", total_size_gb, total_size_gb/(t_end - t_start));
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
    nc_type type;
    int ndims;
    int natts;

    NCERR(nc_open(in_paths[0], NC_NOWRITE, &ncid0));
    NCERR(nc_inq_nvars(ncid0, &nvars));
    
    for (int v=0; v<nvars; ++v) {
	char varname[NC_MAX_NAME+1];
        int ndims;
	NCERR(nc_inq_var(ncid0, v, varname, NULL, &ndims, NULL, NULL));

        int storage0;
        size_t chunk0[ndims];
        int shuffle0;
        int deflate0;
        int deflate_level0;
        NCERR(nc_inq_var_chunking(ncid0, v, &storage0, chunk0));
        NCERR(nc_inq_var_deflate(ncid0, v, &shuffle0, &deflate0, &deflate_level0));

	// fprintf(stdout, "Checking chunking matches for variable \n", varname);

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

		  if (! chunk[d] == chunk0[d]) {
		    fprintf(stderr,"Incompatible chunking for variable %s in %s: %zd - %zd",
			    varname, in_paths[i], chunk[d], chunk0[d]);
		    exit(1);
		  }
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

    if (comm_rank == 0) {
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
            copy_chunked(argv+arg_index+my_file_idx, 1);

            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, current_file_win);
            MPI_Fetch_and_op(&increment, &my_file_idx, MPI_INT, 0, 0, MPI_SUM, current_file_win);
            MPI_Win_unlock(0, current_file_win);
        }
        close_async(0);
    }

    MPI_Win_free(&current_file_win);
    return MPI_Finalize();
}
