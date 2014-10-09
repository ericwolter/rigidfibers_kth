# by default the final executable is next to the source file, to clearly 
# separate the binary from the source lets put it into its own directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

list(APPEND CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wfatal-errors ${OpenMP_C_FLAGS} ${OpenMP_CXX_FLAGS}")

list(APPEND CUDA_NVCC_FLAGS -O3 -gencode arch=compute_30,code=sm_30)

INCLUDE_DIRECTORIES(
    # ${CMAKE_SOURCE_DIR}/lib/local/openblas/include
    # ${CMAKE_SOURCE_DIR}/lib/local/magma/include
    ${CMAKE_SOURCE_DIR}/lib
)

# defines all source code files
set(SOURCES
    parameters.cc
    performance.cc
    resources.cc
    main.cu
    simulation.cu
)

# defines all header code files
set(HEADERS
    common.h
    fiberopt.h
    parameters.h
    resources.h
    simulation.h
)

set(KERNELS
    kernels/assemble_system.cu
    kernels/compute_inner_integral_analytically.cu
    kernels/compute_inner_integral_numerically.cu
    kernels/constants.cu
    kernels/eye_matrix.cu
    kernels/update_fibers.cu
    kernels/update_fibers_firststep.cu
    kernels/update_velocities.cu
)

set(ADDONS
    ../run.py
    ../tools/validate.py
    ../fortran/Int_analytical.f90
    ../fortran/advect_fibers.f90
    ../fortran/coeffs_and_vels.f90
    ../fortran/matrix_assembly.f90
    ../fortran/quadrature.f90
    ../fortran/write_MATLAB.f90
)


# add an executable/binary to the project
CUDA_ADD_EXECUTABLE(fibers ${SOURCES})

# it is not necessary to include the headers here for cmake, however this way
# they also appear in the file browser of different IDEs like Xcode/VisualStudio
add_custom_target(fibers_files DEPENDS ${SOURCES} SOURCES ${HEADERS} ${KERNELS} ${ADDONS})
# CUDA_ADD_CUBLAS_TO_TARGET(fibers)

# TARGET_LINK_LIBRARIES(fibers ${CMAKE_SOURCE_DIR}/lib/local/magma/lib/libmagma.a ${CMAKE_SOURCE_DIR}/lib/local/openblas/lib/libopenblas.a)
