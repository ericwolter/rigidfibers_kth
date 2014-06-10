# by default the final executable is next to the source file, to clearly 
# separate the binary from the source lets put it into its own directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# make compiler as strict as possible, if needed this can be relaxed later
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wpedantic -Wextra -Weverything -Werror -Wfatal-errors")

# defines all source code files
set(SOURCE
    main.cc
    ocl/clplatform.cc
    ocl/cldevice.cc
)

# defines all header code files
set(HEADER
    fiberopt.h
    common.h
    ocl/clplatform.h
    ocl/cldevice.h
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