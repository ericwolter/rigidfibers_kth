#ifndef FIBERS_SIMULATION_H_
#define FIBERS_SIMULATION_H_
/*
 *  simulation.h - header for simulation.cc
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
#include <string>
#include <map>

// Don't bother with warnings in external dependencies
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wdeprecated"
#pragma GCC diagnostic ignored "-Wdocumentation-unknown-command"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wpadded"
#pragma GCC diagnostic ignored "-Wweak-vtables"
#pragma GCC diagnostic ignored "-Wcovered-switch-default"
#pragma GCC diagnostic ignored "-Wextra-semi"
#pragma GCC diagnostic ignored "-Wdocumentation"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#define VIENNACL_WITH_OPENCL
#include "viennacl/ocl/backend.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/lu.hpp"
#pragma GCC diagnostic pop

#include "common.h"
#include "parameters.h"
#include "performance.h"
#include "ocl/cldevice.h"

class Simulation
{
public:
    Simulation(cl_context context, const CLDevice *device, Configuration configuration);
    ~Simulation();

    void step();

    void exportPerformanceMeasurments();
private:
    DISALLOW_COPY_AND_ASSIGN(Simulation);

    cl_context context_;
    const CLDevice* device_;

    Performance* performance_;

    Configuration configuration_;
    size_t global_work_size_;

    cl_command_queue queue_;
    cl_program program_;

    cl_mem previous_position_buffer_;
    cl_mem current_position_buffer_;
    cl_mem next_position_buffer_;
    cl_mem previous_orientation_buffer_;
    cl_mem current_orientation_buffer_;
    cl_mem next_orientation_buffer_;

    cl_mem translational_velocity_buffer_;
    cl_mem rotational_velocity_buffer_;

    cl_mem a_matrix_buffer_;
    cl_mem b_vector_buffer_;
    cl_mem x_vector_buffer_;

    cl_mem quadrature_points_buffer_;
    cl_mem quadrature_weights_buffer_;
    cl_mem legendre_polynomials_buffer_;

    std::map<std::string,cl_kernel> kernels_;

    void initializeQueue();
    void initializeKernels();
    void initializeProgram();
    void initializeBuffers();
    void initializeViennaCL();

    void writeFiberStateToDevice();
    void readFiberStateFromDevice();

    fiberfloat calculateLegendrePolynomial(fiberfloat x, fiberuint n);
    void precomputeLegendrePolynomials();

    void assembleSystem();
    void solveSystem();
    void updateVelocities();

    void dumpLinearSystem();
    void dumpVelocities();
};

#endif // FIBERS_SIMULATION_H_
