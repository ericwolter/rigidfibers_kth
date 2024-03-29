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
#include "parameters.h"
#include "simulation.h"

int main(int argc, char *argv[])
{
    FiberArgs args = fiberopt(argc, argv,/* help */  1, /* version */ "v0.3.0");

    int nDevices;

    cudaGetDeviceCount(&nDevices);
    for (int i = 0; i < nDevices; ++i)
    {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("Device Number: %d\n", i);
        printf("  Device name: %s\n", prop.name);
        printf("  Memory Clock Rate (KHz): %d\n",
               prop.memoryClockRate);
        printf("  Memory Bus Width (bits): %d\n",
               prop.memoryBusWidth);
        printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
               2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6);
    }

    Configuration configuration = Parameters::parseConfigurationFiles(args.parameters, args.layout);

    Parameters::dump(configuration.parameters);

    Simulation simulation(configuration);

    bool running = true;
    unsigned long current_timestep = 0;
    do
    {
        std::cout << "     [CPU]      : Timestep " << current_timestep + 1 << " of " << configuration.parameters.num_timesteps << std::endl;
        simulation.step(current_timestep);

        current_timestep++;

        if (current_timestep >= 1)
        {
            running = false;
        }
    }
    while (running);

    // cleanup
    delete[] configuration.initial_positions;
    delete[] configuration.initial_orientations;
}
