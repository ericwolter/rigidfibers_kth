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
    std::string eventName;
    cl_event event;

    cl_ulong last_time;
    cl_ulong avg_time;
} PerformanceEvent;

class Performance
{
public:
    Performance();
    ~Performance();

    cl_event* getEvent(std::string eventName);
    void updateEvent(std::string eventName);
    void printEvent(std::string eventName);
    void dump();
private:
    DISALLOW_COPY_AND_ASSIGN(Performance);

    std::map<std::string, PerformanceEvent> events_;
};

#endif // FIBERS_PERFORMANCE_H_
