# required a 'new' version of cmake
# this was simply the major.minor version installed by homebrew on the mac at
# the time of the start of the project
cmake_minimum_required(VERSION 2.8)

# the name of this project
# this is simply an internal cmake name to identify this project to for example
# correctly link it into other programs
project(fibers)
enable_language (Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# by default the final executable is next to the source file, to clearly
# separate the binary from the source lets put it into its own directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#list(APPEND CMAKE_Fortran_FLAGS "-cpp -fopenmp -Og -g -fbounds-check -ftree-vectorizer-verbose=2 -ffree-form -ffree-line-length-none")
list(APPEND CMAKE_Fortran_FLAGS "-cpp -fopenmp -O3 -ffast-math -funroll-all-loops -mtune=native -ffree-form -ffree-line-length-none")

# defines all source code files
set(SOURCES
  fibers.f95
)

# defines all header code files
set(HEADERS
  constants.incl
)

# add an executable/binary to the project
add_executable(fibers ${SOURCES} ${HEADERS})
target_link_libraries(fibers ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

set_target_properties(fibers PROPERTIES
                      LINKER_LANGUAGE Fortran)
