# by default the final executable is next to the source file, to clearly 
# separate the binary from the source lets put it into its own directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# make compiler as strict as possible, if needed this can be relaxed later
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wpedantic -Wextra -Weverything -Werror -Wfatal-errors")

# defines all source code files
set(SOURCE
    main.cc
    resources.cc
    simulation.cc
    ocl/cldevice.cc
    ocl/clplatform.cc
    ocl/clutils.cc
)

# defines all header code files
set(HEADER
    common.h
    fiberopt.h
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
)

foreach(KERNEL ${KERNELS})
    get_filename_component(FILENAME ${KERNEL} NAME)
    set(SRC "${CMAKE_SOURCE_DIR}/src/${KERNEL}")
    set(DST "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kernels/${FILENAME}")

    # cmake -E copy_if_different automatically creates the required directories
    # CAUTION:  copy_if_different does not appear to work reliably sometimes the
    #           files are not copied even though they changed.
    #           As a workaround we now ALWAYS copy the files, as long as we
    #           don't have massive amounts of kernels we should be fine...
    add_custom_command(TARGET fibers POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
    )
endforeach(KERNEL)
