/* 
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

#include "stdlib.h"
#include "stdbool.h"

#ifndef READERS_H
#define READERS_H
#ifdef __cplusplus
extern "C" {
#endif

bool get_collated_dim_decomp(int ncid, const char * varname, int decomposition[4]);

//! Get the global length of a collated variable
/**
 *  \param ncid: NetCDF4 file handle
 *  \param varname: Variable name
 *  \returns: Collated variable length
 */
size_t get_collated_dim_len(int ncid, const char * varname);

//! Get collation info from a variable
/**
 * \param ncid: NetCDF4 file handle
 * \param varid: NetCDF4 variable handle
 * \param out_offset[ndims] (out): The offset in the collated array of this file's data
 * \param local_size[ndims] (out): This file's data size
 * \param total_size[ndims] (out): The total collated size of this variable
 * \param ndims: Variable dimensions
 * \returns true if any of the dimensions are collated
 */
bool get_collation_info(int ncid, int varid,
                        size_t out_offset[], size_t local_size[], size_t total_size[],
                        int ndims);

//! Check if any of the dimensions of a variable are collated
/**
 * \param ncid: NetCDF4 file handle
 * \param varid: NetCDF4 variable handle
 * \returns true if any dimension of varid is collated, false otherwise
 */
bool is_collated(int ncid, int varid);

//! Copy all collated variables to the Writer
/**
 *  The main function for Reader processes. Iterates over all collated
 *  variables in the file, choosing to send each variable in compressed (HDF5)
 *  or uncompressed (NetCDF4) mode to the Writer.
 *
 *  \param filename: Input filename
 *  \param async_writer_rank: MPI rank of the writer process
 */
void copy_chunked(const char * filename, int async_writer_rank);

#ifdef __cplusplus
}
#endif
#endif
