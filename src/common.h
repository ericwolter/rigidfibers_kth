#ifndef FIBERS_COMMON_H_
#define FIBERS_COMMON_H_

#ifdef __OPENCL_VERSION__ // device

    // take precautions to easily switch between single/double precision
    // floating points
    #ifdef USE_DOUBLE_PRECISION

        #if defined(cl_khr_fp64)
            #pragma OPENCL EXTENSION cl_khr_fp64 : enable
        #elif defined(cl_amd_fp64)
            #pragma OPENCL EXTENSION cl_amd_fp64 : enable
        #else
            #error "Device does not support double precision extensions!"
        #endif // cl_khr_fp64

        typedef double fiberfloat;
        typedef double4 fiberfloat4;

    #else
        typedef float fiberfloat;
        typedef float4 fiberfloat4;

    #endif // USE_DOUBLE_PRECISION
    
#else // host
    #include <cstdlib>
    #include <iostream>
    
    // Apple uses the name of their platform specifc frameworks as the include
    // path, so we have to manually take care of that
    #if defined(__APPLE__) || defined(__MACOSX)
        #include <OpenCL/opencl.h>
    #else // Windows/Unix
        #include <CL/opencl.h>
    #endif

    #ifdef USE_DOUBLE_PRECISION
        typedef cl_double fiberfloat;
        typedef cl_double4 fiberfloat4;
    #else
        typedef cl_float fiberfloat;
        typedef cl_float4 fiberfloat4;
    #endif // USE_DOUBLE_PRECISION

    // A macro to disallow the copy constructor and operator= functions
    // This should be used in the private: declarations for a class
    // see: https://google-styleguide.googlecode.com/svn/trunk/cppguide.xml?showone=Copy_Constructors#Copy_Constructors
    #define DISALLOW_COPY_AND_ASSIGN(TypeName) \
        TypeName(const TypeName&);             \
        void operator=(const TypeName&)

    inline void clCheckError(cl_int err, const char *name)
    {
        if(err != CL_SUCCESS) 
        {
            std::cerr << "CL ERROR: " << name << " (" << err << ")" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
#endif

#endif // FIBERS_COMMON_H_
