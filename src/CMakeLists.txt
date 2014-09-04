# by default the final executable is next to the source file, to clearly 
# separate the binary from the source lets put it into its own directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# include the libray directory
# 1. ViennaCL is a linear algebra library for many-core architectures
#    It is a header-only library, so we have to simply in include it here
include_directories(
    ${CMAKE_SOURCE_DIR}/lib
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
    # make compiler as strict as possible, if needed this can be relaxed later
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -Wpedantic -Wextra -Weverything -Wfatal-errors")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # using Visual Studio C++
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")
endif()

# defines all source code files
set(SOURCE
    main.cc
    resources.cc
    parameters.cc
    performance.cc
    simulation.cc
    ocl/cldevice.cc
    ocl/clplatform.cc
    ocl/clutils.cc
)

# defines all header code files
set(HEADER
    common.h
    fiberopt.h
    parameters.h
    performance.h
    resources.h
    simulation.h
    ocl/cldevice.h
    ocl/clplatform.h
    ocl/clutils.h
)

# # ensure that all the required header files of the linked libraries can be found
# include_directories(
#     ${OPENCL_INCLUDE_DIRS}
# )

# add an executable/binary to the project
# it is not necessary to include the headers here for cmake, however this way
# they also appear in the file browser of different IDEs like Xcode/VisualStudio
add_executable(fibers ${SOURCE} ${HEADER})

# link correct libraries required for each platform to the executable
if(APPLE)
    target_link_libraries(fibers
        ${OPENCL_LIBRARY}
    )
elseif(MSVC)
    target_link_libraries(fibers
        ${OPENCL_LIBRARY}
    )
else() # Unix
    target_link_libraries(fibers
        ${OPENCL_LIBRARY}
    )
endif()

set(KERNELS
    common.h
    kernels/vadd.cl
    kernels/assemble_system.cl
)

add_custom_target(copy)
foreach(KERNEL ${KERNELS})
    get_filename_component(FILENAME ${KERNEL} NAME)
    set(SRC "${CMAKE_SOURCE_DIR}/src/${KERNEL}")
    set(DST "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kernels/${FILENAME}")

    # cmake -E copy_if_different automatically creates the required directories
    add_custom_command(TARGET copy PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC} ${DST}
    )
endforeach(KERNEL)

add_dependencies(fibers copy)
