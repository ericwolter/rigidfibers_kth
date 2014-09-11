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

#include <stdio.h>

#include <cmath>
#include <ctime>
#include <vector>
#include <sstream>
#include <iomanip>  // used for standard output manipulation (e.g setprecision)

#include "resources.h"

#include "kernels/assemble_system.cu"

Simulation::Simulation(Configuration configuration)
{
    configuration_ = configuration;
    performance_ = new Performance();

    global_work_size_ = IntCeil(configuration_.parameters.num_fibers, 256);

    initializeGPUMemory();

    writeFiberStateToDevice();
    precomputeLegendrePolynomials();
}

Simulation::~Simulation()
{
}

void Simulation::initializeGPUMemory()
{
    const size_t N = configuration_.parameters.num_fibers;
    checkCuda(cudaMalloc(&gpu_previous_positions_, N * sizeof(fiberfloat4)));
    checkCuda(cudaMalloc(&gpu_current_positions_, N * sizeof(fiberfloat4)));
    checkCuda(cudaMalloc(&gpu_next_positions_, N * sizeof(fiberfloat4)));

    checkCuda(cudaMalloc(&gpu_previous_orientations_, N * sizeof(fiberfloat4)));
    checkCuda(cudaMalloc(&gpu_current_orientations_, N * sizeof(fiberfloat4)));
    checkCuda(cudaMalloc(&gpu_next_orientations_, N * sizeof(fiberfloat4)));

    checkCuda(cudaMalloc(&gpu_previous_translational_velocities_, N * sizeof(fiberfloat4)));
    checkCuda(cudaMalloc(&gpu_current_translational_velocities_, N * sizeof(fiberfloat4)));

    checkCuda(cudaMalloc(&gpu_previous_rotational_velocities_, N * sizeof(fiberfloat4)));
    checkCuda(cudaMalloc(&gpu_current_rotational_velocities_, N * sizeof(fiberfloat4)));

    fiberuint num_matrix_rows =
        3 * configuration_.parameters.num_fibers * configuration_.parameters.num_terms_in_force_expansion;
    fiberuint num_matrix_columns = num_matrix_rows;

    checkCuda(cudaMalloc(&gpu_a_matrix_, num_matrix_rows * num_matrix_columns * sizeof(fiberfloat)));
    checkCuda(cudaMalloc(&gpu_b_vector_, num_matrix_rows * sizeof(fiberfloat)));
    // @TODO might be able to just use B vector for the solution
    checkCuda(cudaMalloc(&gpu_x_vector_, num_matrix_rows * sizeof(fiberfloat)));
}

void Simulation::writeFiberStateToDevice()
{
    std::cout << "[CPU] --> [GPU] : Writing initial fiber positions..." << std::endl;
    checkCuda(cudaMemcpy(gpu_current_positions_, configuration_.initial_positions, configuration_.parameters.num_fibers * sizeof(fiberfloat4), cudaMemcpyHostToDevice));
    std::cout << "[CPU] --> [GPU] : Writing initial fiber orientations..." << std::endl;
    checkCuda(cudaMemcpy(gpu_current_orientations_, configuration_.initial_orientations, configuration_.parameters.num_fibers * sizeof(fiberfloat4), cudaMemcpyHostToDevice));
}

void Simulation::readFiberStateFromDevice()
{

}

double Simulation::calculateLegendrePolynomial(double x, unsigned int n)
{
    switch (n)
    {
    case 0:
        return 1;
    case 1:
        return x;
    case 2:
        return (1.0 / 2.0) * (3.0 * pow(x, 2) - 1.0);
    case 3:
        return (1.0 / 2.0) * (5.0 * pow(x, 3) - 3.0 * x);
    case 4:
        return (1.0 / 8.0) * (35.0 * pow(x, 4) - 30.0 * pow(x, 2) + 3.0);
    case 5:
        return (1.0 / 8.0) * (63.0 * pow(x, 5) - 70.0 * pow(x, 3) + 15.0 * x);
    case 6:
        return (1.0 / 16.0) * (231.0 * pow(x, 6) - 315.0 * pow(x, 4) + 105.0 * pow(x, 2) - 5.0);
    case 7:
        return (1.0 / 16.0) * (429.0 * pow(x, 7) - 693.0 * pow(x, 5) + 315.0 * pow(x, 3) - 35.0 * x);
    case 8:
        return (1.0 / 128.0) * (6435.0 * pow(x, 8) - 12012.0 * pow(x, 6) + 6930.0 * pow(x, 4) - 1260.0 * pow(x, 2) + 35.0);
    default:
        std::cerr << "Could not precompute legendre polynomials - n not in range [1..8]: " << n << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Simulation::precomputeLegendrePolynomials()
{
    const size_t total_number_of_points =
        configuration_.parameters.num_quadrature_points_per_interval
        * configuration_.parameters.num_quadrature_intervals;

    // These are the precalculated points for a 3rd order gaussian quadrature
    // These can be looked up in the literature
    double p0 = -sqrt(15.0) / 5.0;
    double p1 = 0.0;
    double p2 = sqrt(15.0) / 5.0;

    // These are the correcponding weights also found in the literature
    double w0 = 5.0 / 9.0;
    double w1 = 8.0 / 9.0;
    double w2 = 5.0 / 9.0;

    // Intialize lower bound of the current integral to -1. At the start of the
    // subinterval iteration this is the lowest bound of the overall integral
    double lower_bound = -1.0;

    // Calculate the size of a single subinterval. The overall integral bounds
    // are [-1, 1] so the range is 2, which can simply be divided by the number
    // of subintervals.
    double interval_size = 2.0 / configuration_.parameters.num_quadrature_intervals;

    fiberfloat *quadrature_points = new fiberfloat[total_number_of_points];
    fiberfloat *quadrature_weights = new fiberfloat[total_number_of_points];
    //  On wikipedia the mapping from [a, b] to [-1, 1] is done with a factor of
    // (b - a) / 2. However in our case b = a + iv, so the factor would simply
    // be iv / 2.
    // Additionally the point as to be shifted by (a + b) / 2, which for us is
    // (a + a + iv) / 2 = (2 * a * iv) / 2.
    // So if we pull out dividing by 2 we arrive at formula below for the point
    // The weight on wikipedia is also scaled by (b - a) / 2, this being iv / 2
    // for us. If we now plug in iv = 2 / NoQI the factor simply becomes
    // 1 / NoQI. So the weights can simply be divided by the number of
    // subintervals as in the formula below
    for (size_t interval_index = 0; interval_index < configuration_.parameters.num_quadrature_intervals; ++interval_index)
    {
        // @TODO potential micro optimizations as p*, w*, interval_size
        //      number_of_quadrature_intervals are constant they could be
        //      calculated outside the loop, however for clarity we leave and
        //      here right now and precomputing polynomials is not performance
        //      critcal anyway
        // @TODO potential memory savings because weights are the same for each
        //      interval
        size_t interval_start_index = interval_index * configuration_.parameters.num_quadrature_points_per_interval;
        quadrature_points[interval_start_index + 0] = (2.0 * lower_bound + interval_size + p0 * interval_size) / 2.0;
        quadrature_points[interval_start_index + 1] = (2.0 * lower_bound + interval_size + p1 * interval_size) / 2.0;
        quadrature_points[interval_start_index + 2] = (2.0 * lower_bound + interval_size + p2 * interval_size) / 2.0;

        quadrature_weights[interval_start_index + 0] = w0 / configuration_.parameters.num_quadrature_intervals;
        quadrature_weights[interval_start_index + 1] = w1 / configuration_.parameters.num_quadrature_intervals;
        quadrature_weights[interval_start_index + 2] = w2 / configuration_.parameters.num_quadrature_intervals;

        // std::cout << quadrature_points[interval_start_index + 0] << std::endl;
        // std::cout << quadrature_points[interval_start_index + 1] << std::endl;
        // std::cout << quadrature_points[interval_start_index + 2] << std::endl;
        // std::cout << std::endl;
        // std::cout << quadrature_weights[interval_start_index + 0] << std::endl;
        // std::cout << quadrature_weights[interval_start_index + 1] << std::endl;
        // std::cout << quadrature_weights[interval_start_index + 2] << std::endl;

        // Advance to next interval by incrementing the lower bound
        lower_bound += interval_size;
    }

    // write quadrature points and weights to device
    checkCuda(cudaMalloc(&gpu_quadrature_points_, total_number_of_points * sizeof(fiberfloat4)));
    checkCuda(cudaMalloc(&gpu_quadrature_weights_, total_number_of_points * sizeof(fiberfloat4)));

    std::cout << "[CPU] --> [GPU] : Writing precomputed quadrature points..." << std::endl;
    checkCuda(cudaMemcpy(gpu_quadrature_points_, quadrature_points, total_number_of_points * sizeof(fiberfloat), cudaMemcpyHostToDevice));
    std::cout << "[CPU] --> [GPU] : Writing precomputed quadrature weights..." << std::endl;
    checkCuda(cudaMemcpy(gpu_quadrature_weights_, quadrature_weights, total_number_of_points * sizeof(fiberfloat), cudaMemcpyHostToDevice));

    // The output matrix contains the legendre polynomials evaluated at each
    // quadrature point. So for each quadrature point we calculate each
    // legendre polynomial up to the number of terms for the force expansion.
    // The results is a matrix where each row represents a point and each column
    // entry represents a legendre polynomial evaluated at that point.
    // The matrix is in column major order as is the default for GLM and GLSL.
    fiberfloat *legendre_polynomials = new fiberfloat[configuration_.parameters.num_terms_in_force_expansion * total_number_of_points];
    for (size_t column_index = 0; column_index < configuration_.parameters.num_terms_in_force_expansion; ++column_index)
    {
        for (size_t point_index = 0; point_index < total_number_of_points; ++point_index)
        {
            legendre_polynomials[point_index + column_index * total_number_of_points] = calculateLegendrePolynomial(quadrature_points[point_index], column_index + 1);
            // std::cout << legendre_polynomials[point_index + column_index * total_number_of_points] << std::endl;
        }
        // std::cout << std::endl;
    }

    // write legendre polynomials to device
    checkCuda(cudaMalloc(&gpu_legendre_polynomials_, sizeof(fiberfloat) * configuration_.parameters.num_terms_in_force_expansion * total_number_of_points));

    std::cout << "[CPU] --> [GPU] : Writing precomputed legendre polynomials..." << std::endl;
    checkCuda(cudaMemcpy(gpu_legendre_polynomials_, legendre_polynomials, sizeof(fiberfloat) * configuration_.parameters.num_terms_in_force_expansion * total_number_of_points, cudaMemcpyHostToDevice));

    // cleanup
    delete[] quadrature_points;
    delete[] quadrature_weights;
    delete[] legendre_polynomials;
}

void Simulation::step(size_t current_timestep)
{
    std::cout << "     [GPU]      : Assembling system..." << std::endl;
    assembleSystem();
    // std::cout << "     [GPU]      : Solving system..." << std::endl;
    // solveSystem();
    // std::cout << "     [GPU]      : Updating velocities..." << std::endl;
    // updateVelocities();

    // //dumpLinearSystem();
    // //dumpVelocities();

    // std::cout << "     [GPU]      : Updating fibers..." << std::endl;
    // updateFibers(current_timestep == 0);

    // DoubleSwap(cl_mem, previous_translational_velocity_buffer_, current_translational_velocity_buffer_);
    // DoubleSwap(cl_mem, previous_rotational_velocity_buffer_, current_rotational_velocity_buffer_);

    // TripleSwap(cl_mem, previous_position_buffer_, current_position_buffer_, next_position_buffer_);
    // TripleSwap(cl_mem, previous_orientation_buffer_, current_orientation_buffer_, next_orientation_buffer_);

    //dumpFibers();
}

void Simulation::assembleSystem()
{
    performance_->start("assemble_system");
    assemble_system <<< (configuration_.parameters.num_fibers + 255) / 256, 256 >>> (
        gpu_current_positions_,
        gpu_current_orientations_,
        gpu_a_matrix_,
        gpu_b_vector_,
        gpu_quadrature_points_,
        gpu_quadrature_weights_,
        gpu_legendre_polynomials_,
        configuration_.parameters.num_fibers
    );
    performance_->stop("assemble_system");
    performance_->print("assemble_system");
}

void Simulation::solveSystem()
{
    // fiberuint num_matrix_rows =
    //     3 * configuration_.parameters.num_fibers * configuration_.parameters.num_terms_in_force_expansion;
    // fiberuint num_matrix_columns = num_matrix_rows;

    // viennacl::matrix<fiberfloat, viennacl::column_major> a_matrix_vienna(a_matrix_buffer_, num_matrix_rows, num_matrix_columns);
    // viennacl::vector<fiberfloat> b_vector_vienna(b_vector_buffer_, num_matrix_rows);
    // viennacl::vector<fiberfloat> x_vector_vienna(x_vector_buffer_, num_matrix_rows);

    // viennacl::linalg::gmres_tag custom_gmres(1e-5, 100, 10);
    // performance_->start("solve_system", true);
    // x_vector_vienna = viennacl::linalg::solve(a_matrix_vienna, b_vector_vienna, custom_gmres);
    // performance_->stop("solve_system");
    // performance_->print("solve_system");
}

void Simulation::updateVelocities()
{
    // cl_int err = 0;

    // cl_uint param = 0; cl_kernel kernel = kernels_["update_velocities"];
    // err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_position_buffer_);
    // err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_orientation_buffer_);
    // err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &x_vector_buffer_);
    // err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_translational_velocity_buffer_);
    // err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_rotational_velocity_buffer_);
    // err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &quadrature_points_buffer_);
    // err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &quadrature_weights_buffer_);
    // err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &legendre_polynomials_buffer_);
    // clCheckError(err, "Could not set kernel arguments for updating velocities");

    // performance_->start("update_velocities", false);
    // err = clEnqueueNDRangeKernel(queue_, kernel, 1, NULL, &global_work_size_, NULL, 0, NULL, performance_->getDeviceEvent("update_velocities"));
    // clCheckError(err, "Could not enqueue kernel");

    // performance_->stop("update_velocities");
    // performance_->print("update_velocities");
}

void Simulation::updateFibers(bool first_timestep)
{
    // // A second order multi-step method
    // // @TODO Why? Which one?
    // // The first time step is a simple forward euler

    // if (first_timestep)
    // {
    //     cl_int err = 0;

    //     cl_uint param = 0; cl_kernel kernel = kernels_["update_fibers_firststep"];
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_position_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &next_position_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_orientation_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &next_orientation_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_translational_velocity_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_rotational_velocity_buffer_);
    //     clCheckError(err, "Could not set kernel arguments for updating fibers");

    //     performance_->start("update_fibers_firststep", false);
    //     err = clEnqueueNDRangeKernel(queue_, kernel, 1, NULL, &global_work_size_, NULL, 0, NULL, performance_->getDeviceEvent("update_fibers_firststep"));
    //     clCheckError(err, "Could not enqueue kernel");

    //     performance_->stop("update_fibers_firststep");
    //     performance_->print("update_fibers_firststep");
    // }
    // else
    // {
    //     cl_int err = 0;

    //     cl_uint param = 0; cl_kernel kernel = kernels_["update_fibers"];
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &previous_position_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_position_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &next_position_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &previous_orientation_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_orientation_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &next_orientation_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &previous_translational_velocity_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_translational_velocity_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &previous_rotational_velocity_buffer_);
    //     err |= clSetKernelArg(kernel, param++, sizeof(cl_mem), &current_rotational_velocity_buffer_);
    //     clCheckError(err, "Could not set kernel arguments for updating fibers");

    //     performance_->start("update_fibers", false);
    //     err = clEnqueueNDRangeKernel(queue_, kernel, 1, NULL, &global_work_size_, NULL, 0, NULL, performance_->getDeviceEvent("update_fibers"));
    //     clCheckError(err, "Could not enqueue kernel");

    //     performance_->stop("update_fibers");
    //     performance_->print("update_fibers");
    // }
}

void Simulation::dumpFibers()
{
    // fiberuint num_rows = 4 * configuration_.parameters.num_fibers;

    // fiberfloat *p = new fiberfloat[num_rows];
    // fiberfloat *o = new fiberfloat[num_rows];

    // cl_int err;
    // err = clEnqueueReadBuffer(queue_, current_position_buffer_, CL_TRUE, 0, sizeof(fiberfloat4) * configuration_.parameters.num_fibers, p, 0, NULL, NULL);
    // clCheckError(err, "Could not read from p");
    // err = clEnqueueReadBuffer(queue_, current_orientation_buffer_, CL_TRUE, 0, sizeof(fiberfloat4) * configuration_.parameters.num_fibers, o, 0, NULL, NULL);
    // clCheckError(err, "Could not read from o");

    // std::string executablePath = Resources::getExecutablePath();

    // std::string p_output_path = executablePath + "/positions.out";
    // std::string o_output_path = executablePath + "/orientations.out";

    // std::ofstream p_output_file;
    // std::ofstream o_output_file;

    // p_output_file.open (p_output_path.c_str());
    // o_output_file.open (o_output_path.c_str());

    // p_output_file << std::fixed << std::setprecision(8);
    // o_output_file << std::fixed << std::setprecision(8);

    // int coordinate = 0;
    // for (fiberuint row_index = 0; row_index < num_rows; ++row_index)
    // {
    //     if (coordinate == 3)
    //     {
    //         coordinate = 0;
    //         continue;
    //     }
    //     else
    //     {
    //         coordinate++;
    //     }

    //     fiberfloat p_value = p[row_index];
    //     fiberfloat o_value = o[row_index];
    //     if (p_value < 0)
    //     {
    //         p_output_file << "     " << p_value;
    //     }
    //     else
    //     {
    //         p_output_file << "      " << p_value;
    //     }
    //     if (o_value < 0)
    //     {
    //         o_output_file << "     " << o_value;
    //     }
    //     else
    //     {
    //         o_output_file << "      " << o_value;
    //     }

    //     p_output_file << std::endl;
    //     o_output_file << std::endl;
    // }
    // p_output_file.close();
    // o_output_file.close();

    // delete[] p;
    // delete[] o;
}

void Simulation::dumpLinearSystem()
{
    // fiberuint num_matrix_rows =
    //     3 * configuration_.parameters.num_fibers * configuration_.parameters.num_terms_in_force_expansion;
    // fiberuint num_matrix_columns = num_matrix_rows;

    // fiberfloat *a_matrix = new fiberfloat[num_matrix_rows * num_matrix_columns];
    // fiberfloat *b_vector = new fiberfloat[num_matrix_rows];
    // fiberfloat *x_vector = new fiberfloat[num_matrix_rows];

    // cl_int err;
    // err = clEnqueueReadBuffer(queue_, a_matrix_buffer_, CL_TRUE, 0, sizeof(fiberfloat) * num_matrix_rows * num_matrix_columns, a_matrix, 0, NULL, NULL);
    // clCheckError(err, "Could not read from a_matrix");
    // err = clEnqueueReadBuffer(queue_, b_vector_buffer_, CL_TRUE, 0, sizeof(fiberfloat) * num_matrix_rows, b_vector, 0, NULL, NULL);
    // clCheckError(err, "Could not read from b_vector");
    // err = clEnqueueReadBuffer(queue_, x_vector_buffer_, CL_TRUE, 0, sizeof(fiberfloat) * num_matrix_rows, x_vector, 0, NULL, NULL);
    // clCheckError(err, "Could not read from x_vector");

    // std::string executablePath = Resources::getExecutablePath();

    // std::string a_matrix_output_path = executablePath + "/a_matrix.out";
    // std::string b_vector_output_path = executablePath + "/b_vector.out";
    // std::string x_vector_output_path = executablePath + "/x_vector.out";

    // std::ofstream a_matrix_output_file;
    // std::ofstream b_vector_output_file;
    // std::ofstream x_vector_output_file;

    // a_matrix_output_file.open (a_matrix_output_path.c_str());
    // b_vector_output_file.open (b_vector_output_path.c_str());
    // x_vector_output_file.open (x_vector_output_path.c_str());

    // a_matrix_output_file << std::fixed << std::setprecision(8);
    // b_vector_output_file << std::fixed << std::setprecision(8);
    // x_vector_output_file << std::fixed << std::setprecision(8);

    // for (fiberuint row_index = 0; row_index < num_matrix_rows; ++row_index)
    // {
    //     for (fiberuint column_index = 0; column_index < num_matrix_columns; ++column_index)
    //     {
    //         fiberfloat value = a_matrix[row_index + column_index * num_matrix_rows];
    //         if (value < 0)
    //         {
    //             a_matrix_output_file << "     " << value;
    //         }
    //         else
    //         {
    //             a_matrix_output_file << "      " << value;
    //         }
    //     }

    //     fiberfloat value;
    //     value = b_vector[row_index];
    //     if (value < 0)
    //     {
    //         b_vector_output_file << "     " << value;
    //     }
    //     else
    //     {
    //         b_vector_output_file << "      " << value;
    //     }
    //     value = x_vector[row_index];
    //     if (value < 0)
    //     {
    //         x_vector_output_file << "     " << value;
    //     }
    //     else
    //     {
    //         x_vector_output_file << "      " << value;
    //     }

    //     a_matrix_output_file << std::endl;
    //     b_vector_output_file << std::endl;
    //     x_vector_output_file << std::endl;
    // }
    // a_matrix_output_file.close();
    // b_vector_output_file.close();
    // x_vector_output_file.close();

    // delete[] a_matrix;
    // delete[] b_vector;
    // delete[] x_vector;
}

void Simulation::dumpVelocities()
{
    // fiberuint num_rows = 4 * configuration_.parameters.num_fibers;

    // fiberfloat *t_vel = new fiberfloat[num_rows];
    // fiberfloat *r_vel = new fiberfloat[num_rows];

    // cl_int err;
    // err = clEnqueueReadBuffer(queue_, current_translational_velocity_buffer_, CL_TRUE, 0, sizeof(fiberfloat4) * configuration_.parameters.num_fibers, t_vel, 0, NULL, NULL);
    // clCheckError(err, "Could not read from t_vel");
    // err = clEnqueueReadBuffer(queue_, current_rotational_velocity_buffer_, CL_TRUE, 0, sizeof(fiberfloat4) * configuration_.parameters.num_fibers, r_vel, 0, NULL, NULL);
    // clCheckError(err, "Could not read from r_vel");

    // std::string executablePath = Resources::getExecutablePath();

    // std::string t_vel_output_path = executablePath + "/t_vel.out";
    // std::string r_vel_output_path = executablePath + "/r_vel.out";

    // std::ofstream t_vel_output_file;
    // std::ofstream r_vel_output_file;

    // t_vel_output_file.open (t_vel_output_path.c_str());
    // r_vel_output_file.open (r_vel_output_path.c_str());

    // t_vel_output_file << std::fixed << std::setprecision(8);
    // r_vel_output_file << std::fixed << std::setprecision(8);

    // int coordinate = 0;
    // for (fiberuint row_index = 0; row_index < num_rows; ++row_index)
    // {
    //     if (coordinate == 3)
    //     {
    //         coordinate = 0;
    //         continue;
    //     }
    //     else
    //     {
    //         coordinate++;
    //     }

    //     fiberfloat t_value = t_vel[row_index];
    //     fiberfloat r_value = r_vel[row_index];
    //     if (t_value < 0)
    //     {
    //         t_vel_output_file << "     " << t_value;
    //     }
    //     else
    //     {
    //         t_vel_output_file << "      " << t_value;
    //     }
    //     if (r_value < 0)
    //     {
    //         r_vel_output_file << "     " << r_value;
    //     }
    //     else
    //     {
    //         r_vel_output_file << "      " << r_value;
    //     }

    //     t_vel_output_file << std::endl;
    //     r_vel_output_file << std::endl;
    // }
    // t_vel_output_file.close();
    // r_vel_output_file.close();

    // delete[] t_vel;
    // delete[] r_vel;
}

void Simulation::exportPerformanceMeasurments()
{
    // performance_->exportMeasurements("performance.out");
}
