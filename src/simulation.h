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

#include "common.h"
#include "parameters.h"
#include "ocl/cldevice.h"

class Simulation
{
public:
    Simulation(cl_context context, const CLDevice *device, Configuration configuration);
    ~Simulation();
private:
    DISALLOW_COPY_AND_ASSIGN(Simulation);

    cl_context context_;
    const CLDevice* device_;

    Configuration configuration_;

    cl_command_queue queue_;
    cl_program program_;

    cl_mem a_buffer_;
    cl_mem b_buffer_;
    cl_mem c_buffer_;

    cl_mem quadrature_points_buffer_;
    cl_mem quadrature_weights_buffer_;
    cl_mem legendre_polynomials_buffer_;

    std::map<std::string,cl_kernel> kernels_;

    void initalizeQueue();
    void initalizeKernels();
    void initalizeProgram();
    void initalizeBuffers();

    void writeFiberStateToDevice();
    void readFiberStateFromDevice();

    fiberfloat calculateLegendrePolynomial(fiberfloat x, fiberuint n);
    void precomputeLegendrePolynomials(fiberuint number_of_quadrature_intervals);
};

#endif // FIBERS_SIMULATION_H_
