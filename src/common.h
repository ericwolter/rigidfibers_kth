#ifndef FIBERS_COMMON_H_
#define FIBERS_COMMON_H_

#ifdef __OPENCL_VERSION__ // device

    // take precautions to easily switch between single/double precision
    // floating points
    #ifdef USE_DOUBLE_PRECISION

        #ifdef cl_khr_fp64
            #pragma OPENCL EXTENSION cl_khr_fp64 : enable
        #elif defined(cl_amd_fp64)
            #pragma OPENCL EXTENSION cl_amd_fp64 : enable
        #else
            #error "Double precision floating point not supported by OpenCL implementation."
        #endif

        typedef double fiberfloat;
        typedef double4 fiberfloat4;

    #else
        typedef float fiberfloat;
        typedef float4 fiberfloat4;

    #endif // USE_DOUBLE_PRECISION

    // see comment in host section for why this is needed
    typedef int fiberint;
    typedef uint fiberuint;
    
#else // host
    #include <cstdlib>
    #include <iostream>
    #include <chrono>
    
    // Apple uses the name of their platform specifc frameworks as the include
    // path, so we have to manually take care of that
    #if defined(__APPLE__)
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

    // If we are using the same struct on both the host and device we have to
    // be careful about the structs alignment (also OpenCL structs can't
    // contain for that reason). To ensure that the struct works everywhere
    // we have to use the cl types on the host (i.e. cl_float, cl_uint etc.)
    // So in order to only define the struct once we define our own types here
    // which are automatically defined to be the correct type when they are used
    // on the host or the device respectively
    typedef cl_int fiberint;
    typedef cl_uint fiberuint;

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

    #define IntCeil(num, divider) ((((num) + (divider) - 1) / (divider)) * (divider))
    #define DoubleSwap(t, a, b) { t tmp = a; a = b; b = tmp; }
    #define TripleSwap(t, a, b, c) {t tmp = a; a = b; b = c; c = tmp; } 

    // This struct is used on both the host and the device that's why it uses only
    // our custom types, see details host sections for more details
    typedef struct
    {
        fiberfloat slenderness;                         // the slenderness parameter
                                                        //   epsilon = a/2L (e.g. 0.01)
        fiberfloat timestep;                            // the timestep size (e.g. 0.1)
        fiberuint num_fibers;                           // the number of fibers
        fiberuint num_timesteps;                        // the number of timesteps
        fiberuint num_terms_in_force_expansion;         // the number of terms used for the
                                                        //   force expansion (e.g. 5)
        fiberuint num_quadrature_intervals;             // the number of intervals the
                                                        //   integral is subdivided into
                                                        //   on each subinterval 3 gaussian
                                                        //   quadrature points are used
                                                        //   (e.g. 8)
        fiberuint num_quadrature_points_per_interval;   // the number of points per
                                                        //   quadrature points (e.g. 3)
        fiberint use_analytical_integration;            // NOT YET USED
        fiberint use_direct_solver;                     // NOT YET USED
    } FiberParams;
    
#endif

#endif // FIBERS_COMMON_H_
