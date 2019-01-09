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
void get_collated_dim_len(int ncid, const char * varname, size_t * len);

// Get collation info from a variable
// out_offset[ndims] - The offset in the collated array of this variable
// total_size[ndims] - The total collated size of this variable
// returns true if any of the dimensions are collated
bool get_collation_info(int ncid, int varid,
                        size_t out_offset[], size_t local_size[], size_t total_size[],
                        int ndims);

// Returns true if any of the dimensions are collated
bool is_collated(int ncid, int varid);

// Copy all chunked variables
void copy_chunked(const char * filename, int async_writer_rank);

#ifdef __cplusplus
}
#endif
#endif
