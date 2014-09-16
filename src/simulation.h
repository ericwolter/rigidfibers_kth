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

#define VIENNACL_WITH_CUDA
#include "viennacl/matrix.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/linalg/bicgstab.hpp"

#include "common.h"
#include "parameters.h"
#include "performance.h"

class Simulation
{
public:
    Simulation(Configuration configuration);
    ~Simulation();

    void step(size_t current_timestep);

    void exportPerformanceMeasurments();
private:
    DISALLOW_COPY_AND_ASSIGN(Simulation);

    Performance* performance_;

    Configuration configuration_;
    size_t global_work_size_;

    fiberfloat4 *gpu_previous_positions_;
    fiberfloat4 *gpu_current_positions_;
    fiberfloat4 *gpu_next_positions_;

    fiberfloat4 *gpu_previous_orientations_;
    fiberfloat4 *gpu_current_orientations_;
    fiberfloat4 *gpu_next_orientations_;

    fiberfloat4 *gpu_previous_translational_velocities_;
    fiberfloat4 *gpu_current_translational_velocities_;

    fiberfloat4 *gpu_previous_rotational_velocities_;
    fiberfloat4 *gpu_current_rotational_velocities_;

    fiberfloat *gpu_a_matrix_;
    fiberfloat *gpu_b_vector_;
    fiberfloat *gpu_x_vector_;

    fiberfloat *gpu_quadrature_points_;
    fiberfloat *gpu_quadrature_weights_;
    fiberfloat *gpu_legendre_polynomials_;

    void initializeGPUMemory();

    void writeFiberStateToDevice();
    void readFiberStateFromDevice();

    double calculateLegendrePolynomial(double x, unsigned int n);
    void precomputeLegendrePolynomials();

    void assembleSystem();
    void solveSystem();
    void updateVelocities();
    void updateFibers(bool first_timestep);

    void dumpFibers();
    void dumpLinearSystem();
    void dumpVelocities();
};

#endif // FIBERS_SIMULATION_H_
