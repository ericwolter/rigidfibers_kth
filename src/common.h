#ifndef FIBERS_COMMON_H_
#define FIBERS_COMMON_H_

#ifdef __OPENCL_VERSION__ // device
    
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

    // A macro to disallow the copy constructor and operator= functions
    // This should be used in the private: declarations for a class
    // see: https://google-styleguide.googlecode.com/svn/trunk/cppguide.xml?showone=Copy_Constructors#Copy_Constructors
    #define DISALLOW_COPY_AND_ASSIGN(TypeName) \
        TypeName(const TypeName&);             \
        void operator=(const TypeName&)

    inline void clCheckError(cl_int err, const char *name) {
        if(err != CL_SUCCESS) {
            std::cerr << "CL ERROR: " << name << " (" << err << ")" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
#endif

#endif // FIBERS_COMMON_H_
