#include "../common.h"

__global__ void assemble_system(
    const fiberfloat4 *positions,
    const fiberfloat4 *orientations,
    fiberfloat *a_matrix,
    fiberfloat *b_vector,
    const fiberfloat *quadrature_points,
    const fiberfloat *quadrature_weights,
    const fiberfloat *legendre_polynomials,
    const fiberuint NUMBER_OF_FIBERS
    )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= NUMBER_OF_FIBERS) return;

    if (i < NUMBER_OF_FIBERS) b_vector[i] = NUMBER_OF_FIBERS * positions[i].x + orientations[i].x;
}
