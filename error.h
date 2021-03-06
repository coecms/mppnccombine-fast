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

#ifndef ERROR_H
#define ERROR_H
#ifdef __cplusplus
extern "C" {
#endif

//! NetCDF error handler
/**
 *  If a NetCDF call has errored reports the error and exits
 *
 *  \param x: The return code of a NetCDF library call
 */
#define NCERR(x) handle_nc_error(x, __FILE__, __LINE__)
void handle_nc_error(int err, const char * file, int line);

//! HDF5 error handler
/**
 *  If a HDF5 call has errored reports the error and exits
 *
 *  \param x: The return code of a HDF5 library call
 */
#define H5ERR(x) handle_h5_error(x, __FILE__, __LINE__)
void handle_h5_error(int err, const char * file, int line);

//! C error handler
/**
 *  If a C library call has errored reports the error and exits
 *
 *  \param x: The return code of a C library call
 *  \param message: Error message
 */
#define CERR(x, message) handle_c_error(x, message, __FILE__, __LINE__)
void handle_c_error(int err, const char * message, const char * file, int line);

#define LOG_DEBUG 4
#define LOG_INFO 3
#define LOG_WARNING 2
#define LOG_ERROR 1

//! Set the output level for log messages
/**
 *  \param level: Messages with a level less than or equal to this will be
 *      output
 */
void set_log_level(int level);

//! Send a message to the log
/**
 *  \param level: Message log level
 *  \param message: Message (printf-like format string)
 *  \param ...: Message arguments
 */
void log_message(int level, const char * message, ...);

#ifdef __cplusplus
}
#endif
#endif
