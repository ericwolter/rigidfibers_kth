/*
 *  simulation.cc - contains all logic required for simulating the fibers
 *
 *  Copyright (C) 2014  Eric Wolter <eric.wolter@gmx.de>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
#include "simulation.h"

#include <vector>

#include "resources.h"

Simulation::Simulation(cl_context context, const CLDevice *device)
{
    context_ = context;
    device_ = device;

    initalizeQueue();
    initalizeProgram();
    initalizeKernels();
    initalizeBuffers();

    cl_int err;
    size_t count = 100;

    cl_uint param = 0; cl_kernel kernel = kernels_["vadd"];
    err  = clSetKernelArg(kernel, param++, sizeof(cl_mem), &a_buffer_);
    err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &b_buffer_);
    err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &c_buffer_);
    err |= clSetKernelArg(kernel, param++, sizeof(cl_uint), &count);
    clCheckError(err, "Could not set kernel arguments");

    const size_t global_work_size = IntCeil(count, 32);
    err = clEnqueueNDRangeKernel(queue_, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
    clCheckError(err, "Could not enqueue kernel");

    fiberfloat *a_data = (fiberfloat *)malloc(sizeof(fiberfloat) * count);
    fiberfloat *b_data = (fiberfloat *)malloc(sizeof(fiberfloat) * count);
    fiberfloat *c_data = (fiberfloat *)malloc(sizeof(fiberfloat) * count);
    err = clEnqueueReadBuffer(queue_, a_buffer_, CL_TRUE, 0, sizeof(fiberfloat) * count, a_data, 0, NULL, NULL);
    clCheckError(err, "Could not read from a buffer");
    err = clEnqueueReadBuffer(queue_, b_buffer_, CL_TRUE, 0, sizeof(fiberfloat) * count, b_data, 0, NULL, NULL);
    clCheckError(err, "Could not read from b buffer");
    err = clEnqueueReadBuffer(queue_, c_buffer_, CL_TRUE, 0, sizeof(fiberfloat) * count, c_data, 0, NULL, NULL);
    clCheckError(err, "Could not read from c buffer");

    fiberfloat tmp;
    for (size_t i = 0; i < count; ++i)
    {
        tmp = a_data[i] + b_data[i];
        tmp -= c_data[i];
        printf("%lu tmp %f h_a %f h_b %f h_c %f \n", i, tmp, a_data[i], b_data[i], c_data[i]);
    }

    free(a_data);
    free(b_data);
    free(c_data);
}

Simulation::~Simulation()
{
    clFinish(queue_);

    clReleaseMemObject(a_buffer_);
    clReleaseMemObject(b_buffer_);
    clReleaseMemObject(c_buffer_);

    clReleaseProgram(program_);
    // @todo release all kernels
    clReleaseCommandQueue(queue_);
    queue_ = NULL;
}

void Simulation::initalizeQueue()
{
    cl_int err;

    queue_ = clCreateCommandQueue(context_, device_->id(), 0, &err);
    clCheckError(err, "Could not create command queue");
}

void Simulation::initalizeProgram()
{
    cl_int err;

    const std::string kernel_filenames[] =
    {
        "common.h",
        "vadd.cl",
        ""
    };

    // load kernel sources
    std::vector<std::string> kernel_sources;
    for (int i = 0; kernel_filenames[i] != ""; i++)
    {
        // Read source from disk
        std::string source = Resources::getKernelSource(kernel_filenames[i]);

        // Load into compile list
        kernel_sources.push_back(source);
    }

    // the OpenCL is a C API and thus only supports const char*. However for
    // convenience we use proper std::string everywhere else and only convert
    // to C-Land at the last moment
    std::vector<const char *> cstr_kernel_sources;
    cstr_kernel_sources.reserve(kernel_sources.size());
    for (size_t i = 0; i < kernel_sources.size(); ++i)
    {
        cstr_kernel_sources.push_back(kernel_sources[i].c_str());
    }

    // -Wshorten-64-to-32 for cstr_kernel_sources.str() is totally fine here
    // We have no other option than to conform to the OpenCL API here and
    // regardless 32bit for the number of kernels should be enough anyway...
    program_ = clCreateProgramWithSource(context_, (cl_uint)cstr_kernel_sources.size(), &cstr_kernel_sources[0], NULL, &err);
    clCheckError(err, "Could not create program from sources");

    cl_int buildError = clBuildProgram(program_, 0, NULL, NULL, NULL, NULL);

    size_t size;
    err = clGetProgramBuildInfo(program_, device_->id(), CL_PROGRAM_BUILD_LOG, 0, NULL, &size);
    clCheckError(err, "Could not get size of build log");
    // only show build log if it actually has something to show. Even if we have
    // no build log the size returned by clGetProgramBuildInfo still contains
    // one char. So we can ignore that one char and only show the build log when
    // we have a size larger than 1.
    if (size > 1)
    {
        char *buildLog = (char *)malloc(sizeof(char) * size);
        err = clGetProgramBuildInfo(program_, device_->id(), CL_PROGRAM_BUILD_LOG, size, buildLog, NULL);
        std::cout << buildLog << std::endl;
        free(buildLog);
    }

    clCheckError(buildError, "Could not build program");
}

void Simulation::initalizeKernels()
{
    cl_int err;

    cl_uint num_kernels;
    err = clCreateKernelsInProgram(program_, 0, NULL, &num_kernels);
    clCheckError(err, "Could not get number of kernels in program");

    cl_kernel *kernels = (cl_kernel *)malloc(sizeof(cl_kernel) * num_kernels);
    err = clCreateKernelsInProgram(program_, num_kernels, kernels, NULL);
    clCheckError(err, "Could not get kernels in program");

    size_t size;
    std::map<std::string, cl_kernel> kernel_map;
    for (cl_uint i = 0; i < num_kernels; ++i)
    {
        cl_kernel kernel = kernels[i];

        err = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &size);
        clCheckError(err, "Could not get length of kernel function name");
        char *cstr_function_name = (char *)malloc(sizeof(char) * size);
        err = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, size, cstr_function_name, NULL);
        clCheckError(err, "Could not get kernel function name");

        std::string function_name = cstr_function_name;
        kernel_map[function_name] = kernel;
    }
    free(kernels);

    kernels_ = kernel_map;
}

void Simulation::initalizeBuffers()
{
    cl_int err;
    size_t count = 100;
    a_buffer_ = clCreateBuffer(context_, CL_MEM_READ_ONLY, sizeof(fiberfloat) * count, NULL, NULL);
    b_buffer_ = clCreateBuffer(context_, CL_MEM_READ_ONLY, sizeof(fiberfloat) * count, NULL, NULL);
    c_buffer_ = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, sizeof(fiberfloat) * count, NULL, NULL);

    fiberfloat *a_data = (fiberfloat *)malloc(sizeof(fiberfloat) * count);
    fiberfloat *b_data = (fiberfloat *)malloc(sizeof(fiberfloat) * count);
    for (size_t i = 0; i < count; ++i)
    {
        a_data[i] = i + 1;
        b_data[i] = count - i;
    }

    err = clEnqueueWriteBuffer(queue_, a_buffer_, CL_TRUE, 0, sizeof(fiberfloat) * count, a_data, 0, NULL, NULL);
    clCheckError(err, "Could not write data to a buffer");
    err = clEnqueueWriteBuffer(queue_, b_buffer_, CL_TRUE, 0, sizeof(fiberfloat) * count, b_data, 0, NULL, NULL);
    clCheckError(err, "Could not write data to b buffer");

    free(a_data);
    free(b_data);
}

