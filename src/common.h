#ifndef FIBERS_COMMON_H_
#define FIBERS_COMMON_H_

#include <cstdlib>
#include <iostream>

// includes CUDA Runtime
#include <cuda_runtime.h>

#include <assert.h>
#define DEBUG

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    std::cerr << "CUDA Runtime Error: " << cudaGetErrorString(result) << std::endl;
    assert(result == cudaSuccess);
  }
#endif
  return result;
}

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
// see: https://google-styleguide.googlecode.com/svn/trunk/cppguide.xml?showone=Copy_Constructors#Copy_Constructors
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);             \
    void operator=(const TypeName&)

#define IntCeil(num, divider) ((((num) + (divider) - 1) / (divider)) * (divider))
#define DoubleSwap(t, a, b) { t tmp = a; a = b; b = tmp; }
#define TripleSwap(t, a, b, c) {t tmp = a; a = b; b = c; c = tmp; }

#ifdef USE_DOUBLE_PRECISION
    typedef double fiberfloat;
    typedef double4 fiberfloat4;
#else
    typedef float fiberfloat;
    typedef float4 fiberfloat4;
#endif // USE_DOUBLE_PRECISION

typedef int fiberint;
typedef uint fiberuint;

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
    fiberint use_analytical_integration;            // option which determines how the
                                                    // the inner integral is evaluated
    fiberint use_direct_solver;                     // NOT YET USED
} FiberParams;
    
#endif // FIBERS_COMMON_H_
