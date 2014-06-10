/*
 *  fibers - simulates slender fibers in a fluid.
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

#include <iostream>

#include "fiberopt.h"

#include "common.h"
#include "parameters.h"
#include "simulation.h"

#include "ocl/clutils.h"

int main(int argc, char *argv[])
{
    FiberArgs args = fiberopt(argc, argv,/* help */  1, /* version */ "v1.0.0-alpha");

    FiberParams params = Parameters::parseParameterFile(args.parameters);
    Parameters::dump(params);

    // std::cout << args.parameters << std::endl;
    // std::cout << args.layout << std::endl;

    const CLPlatform *selected_platform = CLUtils::selectPlatform();
    const CLDevice *selected_device = CLUtils::selectDevice(selected_platform);

    cl_context context = CLUtils::createContext(selected_platform, selected_device);

    Simulation simulation(context, selected_device);

    // cl_command_queue queue = clCreateCommandQueue(context, device_id, 0, &err);

    // cl_program program = clCreateProgramWithSource(context, 1, (const char **) & KernelSource, NULL, &err);

    // std::cout << Resources::getKernelSource("vadd.cl") << std::endl;

    // bool running = true;

    // do
    // {
    //     std::cout << context << std::endl;
    //     std::cout << queue << std::endl;
    // }
    // while (running);

    // cleanup
    clReleaseContext(context);

    return 0;
}
