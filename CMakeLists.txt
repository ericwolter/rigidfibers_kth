# required a 'new' version of cmake
# this was simply the major.minor version installed by homebrew on the mac at
# the time of the start of the project
cmake_minimum_required(VERSION 2.8)

# the name of this project
# this is simply an internal cmake name to identify this project to for example
# correctly link it into other programs
project(fibers)

find_package(CUDA REQUIRED)
find_package(OpenMP REQUIRED)

include(FindCUDA)
include(FindOpenMP)

# add the source directory to allow cmake to compile/link the actual code
add_subdirectory(src)

# WARN: Currently NOT used
# add the lib directory to allow cmake to find all dependencies
#add_subdirectory(lib)