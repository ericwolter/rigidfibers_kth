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

Performance::Performance() {}
Performance::~Performance() {}

cl_event* Performance::getEvent(std::string eventName) 
{
    std::map<std::string, PerformanceEvent>::iterator performance_event = events_.find(eventName);

    // create a new event if this is the first time we encounter it
    if(performance_event == events_.end()) 
    {
        PerformanceEvent event;
        event.eventName = eventName;

        events_[eventName] = event;

        performance_event = events_.find(eventName);
    }

    return &(performance_event->second).event;
}

void Performance::updateEvent(std::string eventName) 
{   
    std::map<std::string, PerformanceEvent>::iterator performance_event = events_.find(eventName);

    clWaitForEvents(1, &performance_event->second.event);

    cl_ulong start, end;
    clGetEventProfilingInfo(performance_event->second.event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
    clGetEventProfilingInfo(performance_event->second.event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);

    performance_event->second.last_time = end - start;
}

void Performance::printEvent(std::string eventName) 
{
    std::map<std::string, PerformanceEvent>::iterator performance_event = events_.find(eventName);
    std::cout << "  [BENCHMARK]   : " << performance_event->second.eventName << " took " << performance_event->second.last_time*1e-09 << " sec" << std::endl;
}

void Performance::dump() 
{
    std::map<std::string, PerformanceEvent>::iterator iter;
    for (iter = events_.begin(); iter != events_.end(); ++iter) {
        PerformanceEvent event = iter->second;
        std::cout << "  [BENCHMARK]   : " << event.eventName << " took " << event.last_time*1e-09 << " sec" << std::endl;
    }
}


