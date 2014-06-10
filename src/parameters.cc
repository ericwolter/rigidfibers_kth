/*
 *  parameters.cc - provides services to parse the parameters and layout files
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
#include "parameters.h"

#include <sstream>

const FiberParams Parameters::parseParameterFile(const std::string parameters_filename)
{
    std::cout << "Parsing parameters file: " << parameters_filename << std::endl;

    std::ifstream parameters_file_stream;
    parameters_file_stream.open(parameters_filename.c_str());
    if (!parameters_file_stream)
    {
        std::cerr << "Could not open parameters file: " << parameters_filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::getline(parameters_file_stream, line);

    // @todo currently only old file format is used
    // // check for version in first line
    // if (line.find("#!version 1.0"))
    // {
    //     // this is a file format version 1.0
    // }
    // else
    // {

    // older file formats didn't specifc a version, so if we don't find it
    // in the first line we assume an old format
    // because the old format might already contain valid parameters in its
    // first line we backtrack to the beginning of the stream and let
    // the specialist file format take care of the file as a whole
    parameters_file_stream.seekg(0, parameters_file_stream.beg);
    return Parameters::parseVersion0ParameterFile(parameters_file_stream);
    // }
}
void Parameters::parseInitialLayoutFile(__unused const std::string layout_filename, __unused fiberfloat4* initialPositions, __unused fiberfloat4* initialOrientation, __unused fiberint *num_of_fibers)
{
}

void Parameters::dump(const FiberParams params) 
{
    std::cout << "**************************************************" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "Number of fibers                   : " << params.num_fibers << std::endl;
    std::cout << "Number of timesteps                : " << params.num_timesteps << std::endl;
    std::cout << "Size of timesteps                  : " << params.timestep << std::endl;
    std::cout << "Slenderness                        : " << params.slenderness << std::endl;
    std::cout << "Number of terms in force expansion : " << params.num_terms_in_force_expansion << std::endl;
    std::cout << "Number of quadrature intervals     : " << params.num_quadrature_intervals << std::endl;
    std::cout << "**************************************************" << std::endl;
}

const FiberParams Parameters::parseVersion0ParameterFile(std::ifstream &parameters_file_stream)
{
    std::cout << "Using version 0 file format" << std::endl;
    FiberParams params;

    // the old file format depends on the ordering of the parameters so have to
    // go through manually it line by line
    std::string line;
    int lineIndex = 0;
    while (std::getline(parameters_file_stream, line))
    {
        std::cout << line << std::endl;

        std::istringstream value(line);
        switch (lineIndex)
        {
        case 0: // Label to run
            break;
        case 1: // Label of indata file
            break;
        case 2: // Number of terms in the force expansion
            value >> params.num_terms_in_force_expansion;
            break;
        case 3: // Size of epsilon
            value >> params.slenderness;
            break;
        case 4: // Timestep
            value >> params.timestep;
            break;
        case 5: // No. of timesteps
            value >> params.num_timesteps;
            break;
        case 6: // Save interval (=1 save every time step)
            break;
        case 7: // Number of timesteps saved in one file
            break;
        case 8: // Number of quad.intervals
            value >> params.num_quadrature_intervals;
            break;
        case 9: // Using analytical integration 1 (0 numerical)
            value >> params.use_analytical_integration;
            break;
        case 10: // Using direct solver (0 iterative then give restart=10, maxiter=100, tol=xx)
            value >> params.use_direct_solver;
            break;
        case 11: // Restart parameter
            break;
        case 12: // Max iter
            break;
        case 13: // Tolerance
            break;
        }
        ++lineIndex;
    }

    parameters_file_stream.close();

    return params;
}
const FiberParams Parameters::parseVersion1ParameterFile(__unused std::ifstream &parameters_file_stream) {
    FiberParams params;
    return params;
}
void Parameters::parseVersion0LayoutFile(__unused std::ifstream &layout_file_stream, __unused fiberfloat4* initialPositions, __unused fiberfloat4 *initialOrientation, __unused fiberint *num_of_fibers)
{
}


// const std::string Resources::getKernelSource(const std::string kernel_filename)
// {
//     // Open file
//     std::string kernel_fullpath = Resources::getPathForKernel(kernel_filename);
//     std::ifstream ifs(kernel_fullpath.c_str());
//     if ( !ifs.is_open() )
//     {
//         std::cerr << "Could not open: " << kernel_fullpath << std::endl;
//         exit(EXIT_FAILURE);
//     }

//     // read content
//     return std::string(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
// }
