#ifndef FIBERS_PERFORMANCE_H_
#define FIBERS_PERFORMANCE_H_
/*
 *  performance.h - header for performance.cc
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
#include <map>
#include <string>
#include "common.h"

typedef struct {
    std::string name;
    
    cl_event event;

    std::chrono::high_resolution_clock::time_point host_start;

    double host_last_time;
    double host_average_time;

    double device_last_time;
    double device_average_time;
} PerformanceTracker;

class Performance
{
public:
    Performance(cl_command_queue queue);
    ~Performance();

    cl_event* getDeviceEvent(std::string name);
    void start(std::string name);
    void stop(std::string name);
    void print(std::string name);
    void dump();
private:
    DISALLOW_COPY_AND_ASSIGN(Performance);

    cl_command_queue queue_;

    std::map<std::string, PerformanceTracker> trackers_;
};

#endif // FIBERS_PERFORMANCE_H_
