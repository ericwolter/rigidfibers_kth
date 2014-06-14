#ifndef FIBERS_PARAMETERS_H_
#define FIBERS_PARAMETERS_H_
/*
 *  parameters.h - header for parameters.cc
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
#include "common.h"
#include <fstream>

// The alignment of the configuration struct does not matter so we can safely
// ignore this warning 
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpadded"
// Simple container for the setup configuration which allows the simulation to
// set itself up correctly
typedef struct {
    FiberParams parameters;
    fiberfloat4 *initial_positions;
    fiberfloat4 *initial_orientations;
} Configuration;
#pragma GCC diagnostic pop

class Parameters
{
public:
    static const Configuration parseConfigurationFiles(const std::string parameters_filename, const std::string layout_filename);
    static const FiberParams parseParameterFile(const std::string parameters_filename);
    static void parseInitialLayoutFile(const std::string layout_filename, fiberfloat4** initialPositions, fiberfloat4** initialOrientations, fiberuint *number_of_fibers);

    static void dump(const FiberParams params);
private:
    static const FiberParams parseVersion1ParameterFile(std::ifstream &parameters_file_stream);
    static const FiberParams parseVersion2ParameterFile(std::ifstream &parameters_file_stream);
    static void parseVersion1LayoutFile(std::ifstream &layout_file_stream, fiberfloat4** initialPositions, fiberfloat4** initialOrientations, fiberuint *number_of_fibers);
};

#endif // FIBERS_PARAMETERS_H_
