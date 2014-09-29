#ifndef FIBERS_EYE_MATRIX_KERNEL_
#define FIBERS_EYE_MATRIX_KERNEL_

#include "../common.h"
#include "constants.cu"

__global__ void eye_matrix(
    fiberfloat *a_matrix
)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= NUMBER_OF_FIBERS) return;

    int x_row_index;
    int y_row_index;
    int z_row_index;

    int x_column_index;
    int y_column_index;
    int z_column_index;

    #pragma unroll
    for (int force_index = 0; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
    {
        x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 0;
        y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 1;
        z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 2;

        x_column_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 0;
        y_column_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 1;
        z_column_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 2;

        a_matrix[x_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = 1;
        a_matrix[y_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = 1;
        a_matrix[z_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = 1;
    }

}

#endif //FIBERS_EYE_MATRIX_KERNEL_