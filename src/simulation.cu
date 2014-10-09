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

// #include "magma.h"

#include "resources.h"

#include "kernels/assemble_system.cu"
#include "kernels/update_velocities.cu"
#include "kernels/update_fibers_firststep.cu"
#include "kernels/update_fibers.cu"
#include "kernels/eye_matrix.cu"

Simulation::Simulation(Configuration configuration)
{
    configuration_ = configuration;
    performance_ = new Performance();

    global_work_size_ = IntCeil(NUMBER_OF_FIBERS, 256);

    initializeGPUMemory();

    //magma_init();
    
    writeFiberStateToDevice();
    precomputeLegendrePolynomials();
}

Simulation::~Simulation()
{
    //magma_finalize();
}

void Simulation::initializeGPUMemory()
{
#ifdef VALIDATE
    checkCuda(cudaMalloc(&gpu_validation_, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * 6 * sizeof(int)));
    checkCuda(cudaMemset(gpu_validation_, 0, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * 6 * sizeof(int)));
#endif //VALIDATE

    checkCuda(cudaMalloc(&gpu_previous_positions_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_current_positions_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_next_positions_, NUMBER_OF_FIBERS * sizeof(float4)));

    checkCuda(cudaMalloc(&gpu_previous_orientations_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_current_orientations_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_next_orientations_, NUMBER_OF_FIBERS * sizeof(float4)));

    checkCuda(cudaMalloc(&gpu_previous_translational_velocities_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_current_translational_velocities_, NUMBER_OF_FIBERS * sizeof(float4)));

    checkCuda(cudaMalloc(&gpu_previous_rotational_velocities_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_current_rotational_velocities_, NUMBER_OF_FIBERS * sizeof(float4)));

    checkCuda(cudaMalloc(&gpu_a_matrix_, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * sizeof(float)));
    checkCuda(cudaMalloc(&gpu_b_vector_, TOTAL_NUMBER_OF_ROWS * sizeof(float)));
    // @TODO might be able to just use B vector for the solution
    checkCuda(cudaMalloc(&gpu_x_vector_, TOTAL_NUMBER_OF_ROWS * sizeof(float)));

    std::cout << "     [GPU]      : Resetting system..." << std::endl;
    checkCuda(cudaMemset(gpu_a_matrix_, 0, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * sizeof(float)));

    performance_->start("eye_matrix");
    eye_matrix <<< (NUMBER_OF_FIBERS + 31) / 32, 32 >>> (gpu_a_matrix_);
    
    performance_->stop("eye_matrix");
    performance_->print("eye_matrix");    

    checkCuda(cudaMemset(gpu_b_vector_, 0, TOTAL_NUMBER_OF_ROWS * sizeof(float)));
    checkCuda(cudaMemset(gpu_x_vector_, 0, TOTAL_NUMBER_OF_ROWS * sizeof(float)));
}

void Simulation::writeFiberStateToDevice()
{
    std::cout << "[CPU] --> [GPU] : Writing initial fiber positions..." << std::endl;
    checkCuda(cudaMemcpy(gpu_current_positions_, configuration_.initial_positions, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyHostToDevice));
    std::cout << "[CPU] --> [GPU] : Writing initial fiber orientations..." << std::endl;
    checkCuda(cudaMemcpy(gpu_current_orientations_, configuration_.initial_orientations, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyHostToDevice));
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
    double interval_size = 2.0 / NUMBER_OF_QUADRATURE_INTERVALS;

    float *host_quadrature_points = new float[TOTAL_NUMBER_OF_QUADRATURE_POINTS];
    float *host_quadrature_weights = new float[TOTAL_NUMBER_OF_QUADRATURE_POINTS];
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
    for (size_t interval_index = 0; interval_index < NUMBER_OF_QUADRATURE_INTERVALS; ++interval_index)
    {
        // @TODO potential micro optimizations as p*, w*, interval_size
        //      number_of_quadrature_intervals are constant they could be
        //      calculated outside the loop, however for clarity we leave and
        //      here right now and precomputing polynomials is not performance
        //      critcal anyway
        // @TODO potential memory savings because weights are the same for each
        //      interval
        size_t interval_start_index = interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL;
        host_quadrature_points[interval_start_index + 0] = (2.0 * lower_bound + interval_size + p0 * interval_size) / 2.0;
        host_quadrature_points[interval_start_index + 1] = (2.0 * lower_bound + interval_size + p1 * interval_size) / 2.0;
        host_quadrature_points[interval_start_index + 2] = (2.0 * lower_bound + interval_size + p2 * interval_size) / 2.0;

        host_quadrature_weights[interval_start_index + 0] = w0 / NUMBER_OF_QUADRATURE_INTERVALS;
        host_quadrature_weights[interval_start_index + 1] = w1 / NUMBER_OF_QUADRATURE_INTERVALS;
        host_quadrature_weights[interval_start_index + 2] = w2 / NUMBER_OF_QUADRATURE_INTERVALS;

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

    std::cout << "[CPU] --> [GPU] : Writing precomputed quadrature points..." << std::endl;
    checkCuda(cudaMemcpyToSymbol(quadrature_points, host_quadrature_points, TOTAL_NUMBER_OF_QUADRATURE_POINTS * sizeof(float)));
    std::cout << "[CPU] --> [GPU] : Writing precomputed quadrature weights..." << std::endl;
    checkCuda(cudaMemcpyToSymbol(quadrature_weights, host_quadrature_weights, TOTAL_NUMBER_OF_QUADRATURE_POINTS * sizeof(float)));

    // The output matrix contains the legendre polynomials evaluated at each
    // quadrature point. So for each quadrature point we calculate each
    // legendre polynomial up to the number of terms for the force expansion.
    // The results is a matrix where each row represents a point and each column
    // entry represents a legendre polynomial evaluated at that point.
    // The matrix is in column major order as is the default for GLM and GLSL.
    float *host_legendre_polynomials = new float[NUMBER_OF_TERMS_IN_FORCE_EXPANSION * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
    for (size_t column_index = 0; column_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++column_index)
    {
        for (size_t point_index = 0; point_index < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++point_index)
        {
            host_legendre_polynomials[point_index + column_index * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = calculateLegendrePolynomial(host_quadrature_points[point_index], column_index + 1);
            // std::cout << legendre_polynomials[point_index + column_index * total_number_of_points] << std::endl;
        }
        // std::cout << std::endl;
    }

    std::cout << "[CPU] --> [GPU] : Writing precomputed legendre polynomials..." << std::endl;
    checkCuda(cudaMemcpyToSymbol(legendre_polynomials, host_legendre_polynomials, NUMBER_OF_TERMS_IN_FORCE_EXPANSION * TOTAL_NUMBER_OF_QUADRATURE_POINTS * sizeof(float)));

    // cleanup
    delete[] host_quadrature_points;
    delete[] host_quadrature_weights;
    delete[] host_legendre_polynomials;

    double *host_double_lambda = new double[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
    double *host_double_eigen = new double[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
    float *host_lambda = new float[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
    float *host_eigen = new float[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];

    double c  = log(SLENDERNESS * SLENDERNESS * M_E);
    double d  = -c;
    double e  = 2.0;
    double cc = 1.0;

    host_double_lambda[0] = 2.0;
    host_double_eigen[0] = ((d - e - cc * host_double_lambda[0]) / 2.0) / (d - cc * host_double_lambda[0]);

    host_lambda[0] = host_double_lambda[0];
    host_eigen[0] = host_double_eigen[0];

    for (size_t force_index = 1; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
    {
        host_double_lambda[force_index] = host_double_lambda[force_index - 1] + 2.0 / (force_index + 1);
        host_double_eigen[force_index] = ((d - e - cc * host_double_lambda[force_index]) / 2.0) / (d - cc * host_double_lambda[force_index]);

        // do all calulcations in double precision but cast to the correct GPU precision
        host_lambda[force_index] = host_double_lambda[force_index];
        host_eigen[force_index] = host_double_eigen[force_index];
    }

    std::cout << "[CPU] --> [GPU] : Writing precomputed lambda values..." << std::endl;
    checkCuda(cudaMemcpyToSymbol(lambda, host_lambda, NUMBER_OF_TERMS_IN_FORCE_EXPANSION * sizeof(float)));

    std::cout << "[CPU] --> [GPU] : Writing precomputed eigen values..." << std::endl;
    checkCuda(cudaMemcpyToSymbol(eigen, host_eigen, NUMBER_OF_TERMS_IN_FORCE_EXPANSION * sizeof(float)));

    delete[] host_double_lambda;
    delete[] host_double_eigen;

    delete[] host_lambda;
    delete[] host_eigen;
}

void Simulation::step(size_t current_timestep)
{
    assembleSystem();
    std::cout << "     [GPU]      : Solving system..." << std::endl;
    solveSystem();
    updateVelocities();

#ifdef VALIDATE
    dumpLinearSystem();
#endif //VALIDATE

    // dumpVelocities();

    std::cout << "     [GPU]      : Updating fibers..." << std::endl;
    updateFibers(current_timestep == 0);

    DoubleSwap(float4*, gpu_previous_translational_velocities_, gpu_current_translational_velocities_);
    DoubleSwap(float4*, gpu_previous_rotational_velocities_, gpu_current_rotational_velocities_);

    TripleSwap(float4*, gpu_previous_positions_, gpu_current_positions_, gpu_next_positions_);
    TripleSwap(float4*, gpu_previous_orientations_, gpu_current_orientations_, gpu_next_orientations_);

    //dumpFibers();
}

void Simulation::assembleSystem()
{
    performance_->start("assemble_system");
#ifdef FORCE_1D
    std::cout << "     [GPU]      : Assembling system 1D..." << std::endl;
    assemble_system <<< (NUMBER_OF_FIBERS + 31) / 32, 32 >>> (
                                                           #ifdef VALIDATE
                                                               gpu_validation_,
                                                           #endif //VALIDATE
        gpu_current_positions_,
        gpu_current_orientations_,
        gpu_a_matrix_,
        gpu_b_vector_
    );
#else
    dim3 block_size;
    block_size.x = 8;
    block_size.y = 8;

    dim3 grid_size;
    grid_size.x = (NUMBER_OF_FIBERS + block_size.x-1) / block_size.x;
    grid_size.y = (NUMBER_OF_FIBERS + block_size.y-1) / block_size.y;

    std::cout << "     [GPU]      : Assembling system 2D..." << std::endl;
    assemble_system <<< grid_size, block_size >>> (
                                                           #ifdef VALIDATE
                                                               gpu_validation_,
                                                           #endif //VALIDATE
        gpu_current_positions_,
        gpu_current_orientations_,
        gpu_a_matrix_,
        gpu_b_vector_
    );
#endif
    performance_->stop("assemble_system");
    performance_->print("assemble_system");
}

void Simulation::solveSystem()
{
    viennacl::matrix_base<float, viennacl::column_major> a_matrix_vienna(gpu_a_matrix_, viennacl::CUDA_MEMORY,
                                    TOTAL_NUMBER_OF_ROWS, 0, 1, TOTAL_NUMBER_OF_ROWS,
                                    TOTAL_NUMBER_OF_ROWS, 0, 1, TOTAL_NUMBER_OF_ROWS);
    viennacl::vector<float> b_vector_vienna(gpu_b_vector_, viennacl::CUDA_MEMORY, TOTAL_NUMBER_OF_ROWS);
    viennacl::vector<float> x_vector_vienna(gpu_x_vector_, viennacl::CUDA_MEMORY, TOTAL_NUMBER_OF_ROWS);

    viennacl::linalg::gmres_tag custom_gmres(1e-5, 100, 10);
    performance_->start("solve_system");

    x_vector_vienna = viennacl::linalg::solve(a_matrix_vienna, b_vector_vienna, custom_gmres);

    performance_->stop("solve_system");
    performance_->print("solve_system");

    // magma_int_t *ipiv=NULL;
    // magma_int_t ldda = ((num_matrix_rows+31)/32)*32;  // round up to multiple of 32 for best GPU performance
    // magma_int_t lddx = ldda;
    // magma_int_t info = 0;

    // magma_imalloc_cpu( &ipiv, num_matrix_rows );  // ipiv always on CPU
}

void Simulation::updateVelocities()
{
    performance_->start("update_velocities");
#ifdef FORCE_1D
    std::cout << "     [GPU]      : Updating velocities 1D..." << std::endl;
    update_velocities <<< (NUMBER_OF_FIBERS + 31) / 32, 32 >>> (
        gpu_current_positions_,
        gpu_current_orientations_,
        gpu_x_vector_,
        gpu_current_translational_velocities_,
        gpu_current_rotational_velocities_
    );
#else
    dim3 block_size;
    block_size.x = 8;
    block_size.y = 8;

    dim3 grid_size;
    grid_size.x = (NUMBER_OF_FIBERS + block_size.x-1) / block_size.x;
    grid_size.y = (NUMBER_OF_FIBERS + block_size.y-1) / block_size.y;

    std::cout << "     [GPU]      : Updating velocities 2D..." << std::endl;
    update_velocities <<< grid_size, block_size >>> (
        gpu_current_positions_,
        gpu_current_orientations_,
        gpu_x_vector_,
        gpu_current_translational_velocities_,
        gpu_current_rotational_velocities_
    );
#endif
    performance_->stop("update_velocities");
    performance_->print("update_velocities");    
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
    // A second order multi-step method
    // @TODO Why? Which one?
    // The first time step is a simple forward euler

    if (first_timestep)
    {
        performance_->start("update_fibers_firststep");
        update_fibers_firststep <<< (NUMBER_OF_FIBERS + 31) / 32, 32 >>> (
            gpu_current_positions_,
            gpu_next_positions_,
            gpu_current_orientations_,
            gpu_next_orientations_,
            gpu_current_translational_velocities_,
            gpu_current_rotational_velocities_
        );
        performance_->stop("update_fibers_firststep");
        performance_->print("update_fibers_firststep");
    }
    else
    {
        performance_->start("update_fibers");
        update_fibers <<< (NUMBER_OF_FIBERS + 31) / 32, 32 >>> (
            gpu_previous_positions_,
            gpu_current_positions_,
            gpu_next_positions_,
            gpu_previous_orientations_,
            gpu_current_orientations_,
            gpu_next_orientations_,
            gpu_previous_translational_velocities_,
            gpu_current_translational_velocities_,
            gpu_previous_rotational_velocities_,
            gpu_current_rotational_velocities_
        );
        performance_->stop("update_fibers");
        performance_->print("update_fibers");
    }
}

void Simulation::dumpFibers()
{
    float4 *p = new float4[NUMBER_OF_FIBERS];
    float4 *o = new float4[NUMBER_OF_FIBERS];

    checkCuda(cudaMemcpy(p, gpu_current_positions_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(o, gpu_current_orientations_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));

    std::string executablePath = Resources::getExecutablePath();

    std::string p_output_path = executablePath + "/positions.out";
    std::string o_output_path = executablePath + "/orientations.out";

    std::ofstream p_output_file;
    std::ofstream o_output_file;

    p_output_file.open (p_output_path.c_str());
    o_output_file.open (o_output_path.c_str());

    p_output_file << std::fixed << std::setprecision(8);
    o_output_file << std::fixed << std::setprecision(8);

    for (size_t row_index = 0; row_index < NUMBER_OF_FIBERS; ++row_index)
    {
        float4 p_value = p[row_index];
        float4 o_value = o[row_index];

        p_output_file << (p_value.x < 0 ? "     " : "      ") << p_value.x << std::endl;
        p_output_file << (p_value.y < 0 ? "     " : "      ") << p_value.y << std::endl;
        p_output_file << (p_value.z < 0 ? "     " : "      ") << p_value.z << std::endl;

        o_output_file << (o_value.x < 0 ? "     " : "      ") << o_value.x << std::endl;
        o_output_file << (o_value.y < 0 ? "     " : "      ") << o_value.y << std::endl;
        o_output_file << (o_value.z < 0 ? "     " : "      ") << o_value.z << std::endl;
    }
    p_output_file.close();
    o_output_file.close();

    delete[] p;
    delete[] o;
}

void Simulation::dumpLinearSystem()
{
    float *a_matrix = new float[TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS];
    float *b_vector = new float[TOTAL_NUMBER_OF_ROWS];
    float *x_vector = new float[TOTAL_NUMBER_OF_ROWS];

    checkCuda(cudaMemcpy(a_matrix, gpu_a_matrix_, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * sizeof(float), cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(b_vector, gpu_b_vector_, TOTAL_NUMBER_OF_ROWS * sizeof(float), cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(x_vector, gpu_x_vector_, TOTAL_NUMBER_OF_ROWS * sizeof(float), cudaMemcpyDeviceToHost));

    std::string executablePath = Resources::getExecutablePath();

    std::string a_matrix_output_path = executablePath + "/a_matrix.out";
    std::string b_vector_output_path = executablePath + "/b_vector.out";
    std::string x_vector_output_path = executablePath + "/x_vector.out";

    std::ofstream a_matrix_output_file;
    std::ofstream b_vector_output_file;
    std::ofstream x_vector_output_file;

    a_matrix_output_file.open (a_matrix_output_path.c_str());
    b_vector_output_file.open (b_vector_output_path.c_str());
    x_vector_output_file.open (x_vector_output_path.c_str());

    a_matrix_output_file << std::fixed << std::setprecision(8);
    b_vector_output_file << std::fixed << std::setprecision(8);
    x_vector_output_file << std::fixed << std::setprecision(8);

    for (int row_index = 0; row_index < TOTAL_NUMBER_OF_ROWS; ++row_index)
    {
        for (int column_index = 0; column_index < TOTAL_NUMBER_OF_ROWS; ++column_index)
        {
            float value = a_matrix[row_index + column_index * TOTAL_NUMBER_OF_ROWS];
            if (value < 0)
            {
                a_matrix_output_file << "     " << value;
            }
            else
            {
                a_matrix_output_file << "      " << value;
            }
        }

        float value;
        value = b_vector[row_index];
        if (value < 0)
        {
            b_vector_output_file << "     " << value;
        }
        else
        {
            b_vector_output_file << "      " << value;
        }
        value = x_vector[row_index];
        if (value < 0)
        {
            x_vector_output_file << "     " << value;
        }
        else
        {
            x_vector_output_file << "      " << value;
        }

        a_matrix_output_file << std::endl;
        b_vector_output_file << std::endl;
        x_vector_output_file << std::endl;
    }
    a_matrix_output_file.close();
    b_vector_output_file.close();
    x_vector_output_file.close();

    delete[] a_matrix;
    delete[] b_vector;
    delete[] x_vector;

#ifdef VALIDATE
    int *validation = new int[TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * 6];
    checkCuda(cudaMemcpy(validation, gpu_validation_, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * 6 * sizeof(int), cudaMemcpyDeviceToHost));

    std::string mapping_output_path = executablePath + "/current.map";
    std::ofstream mapping_output_file;

    mapping_output_file.open (mapping_output_path.c_str());

    for (int row_index = 0; row_index < TOTAL_NUMBER_OF_ROWS; ++row_index)
    {
        for (int column_index = 0; column_index < TOTAL_NUMBER_OF_ROWS; ++column_index)
        {
            mapping_output_file << "[V]"
                                << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0]
                                << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1]
                                << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2]
                                << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3]
                                << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4]
                                << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5]
                                << "|" << row_index
                                << "|" << column_index
                                << std::endl;
        }
    }
    mapping_output_file.close();
    delete[] validation;
#endif //VALIDATE
}

void Simulation::dumpVelocities()
{
    float4 *t_vel = new float4[NUMBER_OF_FIBERS];
    float4 *r_vel = new float4[NUMBER_OF_FIBERS];

    checkCuda(cudaMemcpy(t_vel, gpu_current_translational_velocities_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(r_vel, gpu_current_rotational_velocities_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));

    std::string executablePath = Resources::getExecutablePath();

    std::string t_vel_output_path = executablePath + "/t_vel.out";
    std::string r_vel_output_path = executablePath + "/r_vel.out";

    std::ofstream t_vel_output_file;
    std::ofstream r_vel_output_file;

    t_vel_output_file.open (t_vel_output_path.c_str());
    r_vel_output_file.open (r_vel_output_path.c_str());

    t_vel_output_file << std::fixed << std::setprecision(8);
    r_vel_output_file << std::fixed << std::setprecision(8);

    for (size_t row_index = 0; row_index < NUMBER_OF_FIBERS; ++row_index)
    {
        float4 t_value = t_vel[row_index];
        float4 r_value = r_vel[row_index];

        t_vel_output_file << (t_value.x < 0 ? "     " : "      ") << t_value.x << std::endl;
        t_vel_output_file << (t_value.y < 0 ? "     " : "      ") << t_value.y << std::endl;
        t_vel_output_file << (t_value.z < 0 ? "     " : "      ") << t_value.z << std::endl;

        r_vel_output_file << (r_value.x < 0 ? "     " : "      ") << r_value.x << std::endl;
        r_vel_output_file << (r_value.y < 0 ? "     " : "      ") << r_value.y << std::endl;
        r_vel_output_file << (r_value.z < 0 ? "     " : "      ") << r_value.z << std::endl;
    }
    t_vel_output_file.close();
    r_vel_output_file.close();

    delete[] t_vel;
    delete[] r_vel;
}

void Simulation::exportPerformanceMeasurments()
{
    // performance_->exportMeasurements("performance.out");
}
