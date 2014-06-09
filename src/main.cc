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

int main(int argc, char* argv[])
{
    FiberArgs args = fiberopt(argc, argv, /* help */ 1, /* version */ "v1.0.0-alpha");

    printf("Commands\n");
    printf("Arguments\n");
    printf("    input == %s\n", args.input);
    printf("Flags\n");
    printf("    --gui == %s\n", args.gui ? "true" : "false");
    printf("    --help == %s\n", args.help ? "true" : "false");
    printf("    --version == %s\n", args.version ? "true" : "false");
    printf("Options\n");
    return 0;
}
