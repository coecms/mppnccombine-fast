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

#ifndef ASYNC_H
#define ASYNC_H
#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include "hdf5.h"
#include "mpi.h"
#include "netcdf.h"

typedef struct {
    int idx; 
} varid_t;


//! Get info about a variable in the output file
/**
 * For direct chunk writes to work the compression parameters must match
 * between the input and output files. Use this function to query the
 * parameters of the output file, then compare the results against the values
 * returned by nc_inq_var_deflate() on the input file.
 *
 * \param varid: Variable handle obtained with open_variable_async
 * \param ndims: Number of dimensions
 * \param chunk[ndims] (out): Chunk shape
 * \param deflate (out): Compression enabled
 * \param deflate_level (out): Compression level
 * \param shuffle (out): Shuffle filter enabled
 * \param async_writer_rank: MPI rank of writer process
 */
void variable_info_async(
    varid_t varid,
    size_t ndims,
    size_t chunk[],
    int * deflate,
    int * deflate_level,
    int * shuffle,
    int async_writer_rank
    );


//! Write data to the file, using the dataset filters
/**
 * This writes uncompressed data, as obtained by nc_get_vara() on the input
 * file. It is slower than direct chunk writes, as the data must be put through
 * a compression filter, but is more flexible as it can be used to write
 * partial or unaligned chunks.
 *
 * \param varid: Variable handle obtained with open_variable_async
 * \param ndims: Number of dimensions
 * \param offset[ndims]: Offset of this data's origin in the collated dataset
 * \param shape[ndims]: Shape of this data array
 * \param buffer: Compressed chunk data
 * \param type: NetCDF type of the data
 * \param async_writer_rank: MPI rank of writer process
 * \param request (out): MPI request for the communication
 *
 * request must be sent to MPI_Wait for the message to complete
 */
void write_uncompressed_async(
    varid_t varid,
    size_t ndims,
    const size_t offset[],
    const size_t shape[],
    const void * buffer,
    nc_type type,
    int async_writer_rank,
    MPI_Request * request);

//! Write a compressed chunk directly to the file
/**
 * This writes compressed chunks, as obtained by opening the input file in HDF5
 * and calling H5DOread_chunk(). This is faster than copying uncompressed data,
 * but the chunking and compression parameters must be identical on the input
 * and output files, and the chunk must lay on the chunk boundary of the output
 * file. variable_info_async() can be used to determine the chunk layout and
 * compression settings of the variable in the output file.
 *
 * \param varid: Variable handle obtained with open_variable_async
 * \param ndims: Number of dimensions
 * \param filter_mask: HDF5 filter information (must match the output file)
 * \param offset[ndims]: Offset of this chunk's origin in the collated dataset (must
 *             be on a chunk boundary of the output file)
 * \param data_size: Size of the compressed chunk in bytes
 * \param buffer: Compressed chunk data
 * \param async_writer_rank: MPI rank of writer process
 * \param request [out]: MPI request for the communication
 *
 * request must be sent to MPI_Wait for the message to complete
 */
void write_chunk_async(
    varid_t varid,
    size_t ndims,
    uint32_t filter_mask,
    const hsize_t offset[],
    size_t data_size,
    const void * buffer,
    int async_writer_rank,
    MPI_Request * request
    );

//! Open a variable in the async writer
/**
 * \param varname: Variable name
 * \param len: Length of varname, including the closing '/0'
 * \param async_writer_rank: MPI rank of writer process
 *
 * \returns a handle to the variable in the output file
 */
varid_t open_variable_async(
    const char * varname,
    size_t len, // Total length including '/0'
    int async_writer_rank
    );

//! Close a variable in the async writer
/**
 * \param varid: Variable handle obtained with open_variable_async
 * \param async_writer_rank: MPI rank of writer process
 *
 * varid is no longer a valid handle after this call
 */
void close_variable_async(
    varid_t varid,
    int async_writer_rank
    );

//! Close the async writer
/**
 * \param async_writer_rank: MPI rank of writer process
 */
void close_async(
    int async_writer_rank
    );

//! Async runner to accept writes
/**
 * Called by the Writer to accept async messages sent by the Readers. Once all
 * Readers have sent close_async() messages this will return
 *
 * \param filename: Output filename
 *
 * \returns total size written
 */
size_t run_async_writer(
    const char * filename
    );

#ifdef __cplusplus
}
#endif
#endif
