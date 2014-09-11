# by default the final executable is next to the source file, to clearly 
# separate the binary from the source lets put it into its own directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

list(APPEND CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wfatal-errors")

# defines all source code files
set(SOURCES
    parameters.cc
    performance.cc
    main.cu
    simulation.cu
)

# defines all header code files
set(HEADERS
    common.h
    fiberopt.h
    parameters.h
)

# add an executable/binary to the project
# it is not necessary to include the headers here for cmake, however this way
# they also appear in the file browser of different IDEs like Xcode/VisualStudio
CUDA_ADD_EXECUTABLE(fibers ${SOURCES} ${HEADERS})
