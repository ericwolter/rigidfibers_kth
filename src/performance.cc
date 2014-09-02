/*
 *  performance.cc - allows to monitor opencl performance
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
#include "performance.h"

Performance::Performance(cl_command_queue queue)
{
    queue_ = queue;
}
Performance::~Performance() {}

cl_event* Performance::getDeviceEvent(std::string name) 
{
    std::map<std::string, PerformanceTracker>::iterator performance_tracker = trackers_.find(name);

    // create a new event if this is the first time we encounter it
    if(performance_tracker == trackers_.end()) 
    {
        PerformanceTracker tracker;
        tracker.name = name;

        trackers_[name] = tracker;

        performance_tracker = trackers_.find(name);
    }

    return &(performance_tracker->second).event;
}

void Performance::start(std::string name)
{
    clFinish(queue_);
    std::map<std::string, PerformanceTracker>::iterator performance_tracker = trackers_.find(name);

    // create a new event if this is the first time we encounter it
    if(performance_tracker == trackers_.end()) 
    {
        PerformanceTracker tracker;
        tracker.name = name;

        trackers_[name] = tracker;

        performance_tracker = trackers_.find(name);
    }

    performance_tracker->second.host_start = std::chrono::high_resolution_clock::now();
}

void Performance::stop(std::string name)
{
    clFinish(queue_);
    std::chrono::high_resolution_clock::time_point host_end = std::chrono::high_resolution_clock::now();

    std::map<std::string, PerformanceTracker>::iterator performance_tracker = trackers_.find(name);

    std::chrono::high_resolution_clock::time_point host_start = performance_tracker->second.host_start;
    performance_tracker->second.host_last_time = (std::chrono::duration_cast<std::chrono::duration<double> >(host_end - host_start)).count();

    cl_ulong device_start, device_end;
    clGetEventProfilingInfo(performance_tracker->second.event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &device_start, NULL);
    clGetEventProfilingInfo(performance_tracker->second.event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &device_end, NULL);

    performance_tracker->second.device_last_time = (device_end - device_start) * 1e-9;
}

void Performance::print(std::string name) 
{
    std::map<std::string, PerformanceTracker>::iterator performance_tracker = trackers_.find(name);

    std::cout << "  [BENCHMARK]   : " << performance_tracker->second.name << "(host)   took " << performance_tracker->second.host_last_time << " sec" << std::endl;
    std::cout << "  [BENCHMARK]   : " << performance_tracker->second.name << "(device) took " << performance_tracker->second.device_last_time << " sec" << std::endl;
}

void Performance::dump() 
{
    std::map<std::string, PerformanceTracker>::iterator iter;
    for (iter = trackers_.begin(); iter != trackers_.end(); ++iter) {
        print(iter->second.name);
    }
}


