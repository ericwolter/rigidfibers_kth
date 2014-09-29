#ifndef FIBERS_CONSTANTS_KERNEL_
#define FIBERS_CONSTANTS_KERNEL_

#define DIMENSIONS (3)
#define NUMBER_OF_FIBERS (100)
#define TIMESTEP (0.1f)
#define SLENDERNESS (0.01f)
#define NUMBER_OF_TERMS_IN_FORCE_EXPANSION (5)
#define NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL (3)
#define NUMBER_OF_QUADRATURE_INTERVALS (8)
#define TOTAL_NUMBER_OF_QUADRATURE_POINTS (NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL * NUMBER_OF_QUADRATURE_INTERVALS)
// #define USE_ANALYTICAL_INTEGRATION (true)

#define TOTAL_NUMBER_OF_ROWS (NUMBER_OF_FIBERS * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS)

__constant__ fiberfloat quadrature_points[TOTAL_NUMBER_OF_QUADRATURE_POINTS];
__constant__ fiberfloat quadrature_weights[TOTAL_NUMBER_OF_QUADRATURE_POINTS];
__constant__ fiberfloat legendre_polynomials[TOTAL_NUMBER_OF_QUADRATURE_POINTS * NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
__constant__ fiberfloat lambda[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
__constant__ fiberfloat eigen[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];

#endif //FIBERS_CONSTANTS_KERNEL_