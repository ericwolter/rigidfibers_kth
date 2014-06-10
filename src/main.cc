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

#include "common.h"
#include "fiberopt.h"

#include "ocl/clplatform.h"
 
int main(int argc, char* argv[])
{
    FiberArgs args = fiberopt(argc, argv, /* help */ 1, /* version */ "v1.0.0-alpha");

    if(args.gui) {
        
    }

    std::vector<CLPlatform*> platforms = CLPlatform::list();

    if (platforms.size() <= 0) {
        std::cerr << "No platforms found" << std::endl;
        exit(EXIT_FAILURE);
    }

    cl_uint index;

    std::cout << "[?] What platform would you like to use?" << std::endl;
    std::vector<CLPlatform*>::const_iterator itPlatform;
    for (index = 0, itPlatform = platforms.begin(); itPlatform != platforms.end(); ++itPlatform, ++index) {
        CLPlatform *platform = *itPlatform;
        const char* platformName = platform->name();
        const char* platformVendor = platform->vendor();
        std::cout << index << ": " 
                  << platformName << " (" << platformVendor <<")" 
                  << std::endl;
        free((char*)platformName);
        free((char*)platformVendor);
    }

    cl_uint selectedPlatformIndex;

    std::cout << "> ";
    std::cin >> selectedPlatformIndex;

    std::cout << selectedPlatformIndex << std::endl;

    CLPlatform *selectedPlatform = platforms.at(selectedPlatformIndex);
    std::vector<CLDevice*> devices = selectedPlatform->devices();

    if(devices.size() <= 0) {
        std::cerr << "No devices found" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "[?] Which device would you like to use?" << std::endl;
    std::vector<CLDevice*>::const_iterator itDevice;
    for (index = 0, itDevice = devices.begin(); itDevice != devices.end(); ++itDevice, ++index) {
        CLDevice *device = *itDevice;
        const char* deviceName = device->name();
        std::cout << index << ": " 
                  << deviceName
                  << std::endl;
        free((char*)deviceName);
    }

    cl_uint selectedDeviceIndex;

    std::cout << "> ";
    std::cin >> selectedDeviceIndex;

    std::cout << selectedDeviceIndex << std::endl;

    return 0;
}
