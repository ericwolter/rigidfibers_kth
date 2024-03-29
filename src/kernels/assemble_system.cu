#ifndef FIBERS_ASSEMBLE_SYSTEM_KERNEL_
#define FIBERS_ASSEMBLE_SYSTEM_KERNEL_

#include "constants.cu"
#include "compute_inner_integral_analytically.cu"
#include "compute_inner_integral_numerically.cu"

__global__ void assemble_system(
    const fiberfloat4 *positions,
    const fiberfloat4 *orientations,
    fiberfloat *a_matrix,
    fiberfloat *b_vector
)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= NUMBER_OF_FIBERS) return;
    if (j >= NUMBER_OF_FIBERS) return;

    const fiberfloat c  = logf(SLENDERNESS * SLENDERNESS * M_E);
    const fiberfloat d  = -c;
    const fiberfloat e  = 2.0f;
    const fiberfloat cc = 1.0f;
    const fiberfloat D1 = 0.75f / (d - 2.0f * cc);

    const fiberfloat4 position_i = positions[i];
    const fiberfloat4 orientation_i = orientations[i];

    fiberfloat4 external_force;
    external_force.x = 0.5f * 0.0f;
    external_force.y = 0.5f * 0.0f;
    external_force.z = 0.5f * -1.0f;

    fiberuint x_row_index;
    fiberuint y_row_index;
    fiberuint z_row_index;

    fiberuint x_column_index;
    fiberuint y_column_index;
    fiberuint z_column_index;

    // b_vector[i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + 0] = 0.0f;
    // b_vector[i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + 1] = 0.0f;
    // b_vector[i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + 2] = 0.0f;

    fiberfloat G[TOTAL_NUMBER_OF_QUADRATURE_POINTS * 6];
    fiberfloat GF[TOTAL_NUMBER_OF_QUADRATURE_POINTS * 3];

    // for (int force_index = 1; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
    // {
    //     x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 0;
    //     y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 1;
    //     z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 2;

    //     b_vector[x_row_index] = 0.0f;
    //     b_vector[y_row_index] = 0.0f;
    //     b_vector[z_row_index] = 0.0f;
    // }

    if (i == j) return;

    const fiberfloat4 position_j = positions[j];
    const fiberfloat4 orientation_j = orientations[j];

    for (fiberuint force_index_i = 0; force_index_i < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index_i)
    {
        fiberfloat Q1;
        fiberfloat Q2;
        fiberfloat Q3;

        fiberuint force_index_j = 0;

        // theta in equation 23
        fiberfloat T11 = 0.0f;
        fiberfloat T22 = 0.0f;
        fiberfloat T33 = 0.0f;
        fiberfloat T12 = 0.0f;
        fiberfloat T13 = 0.0f;
        fiberfloat T23 = 0.0f;

        fiberfloat TF1 = 0.0f;
        fiberfloat TF2 = 0.0f;
        fiberfloat TF3 = 0.0f;
        fiberfloat QF;

#ifdef USE_ANALYTICAL_INTEGRATION
        compute_G_analytic(position_i, orientation_i, position_j, orientation_j, force_index_i, external_force, G, GF, i == 89 && j == 21);
#else
        compute_G_numeric(position_i, orientation_i, position_j, orientation_j, force_index_i, external_force, G, GF, i == 89 && j == 21);
#endif

        for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
        {
            const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_i];
            const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
            T11 += quadrature_weight * G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            T22 += quadrature_weight * G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            T33 += quadrature_weight * G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            T12 += quadrature_weight * G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            T13 += quadrature_weight * G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            T23 += quadrature_weight * G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;

            if (force_index_i == 0)
            {
                TF1 += quadrature_weight * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                TF2 += quadrature_weight * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                TF3 += quadrature_weight * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            }
        }

        Q1 = T11 * orientation_i.x + T12 * orientation_i.y + T13 * orientation_i.z;
        Q2 = T12 * orientation_i.x + T22 * orientation_i.y + T23 * orientation_i.z;
        Q3 = T13 * orientation_i.x + T23 * orientation_i.y + T33 * orientation_i.z;

        //if (i == 0 && j == 1 && force_index_i == 0)
        //{
        //    printf("i=%d;j=%d;force_index_i=%d\nTF1=%f;TF2=%f;TF3=%f;QF=%f\n", i, j, force_index_i, TF1, TF2, TF3, QF);
        //}

        x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 0;
        y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 1;
        z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 2;

        x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 0;
        y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 1;
        z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 2;

        // if(i==0) {
        //     printf("%d,%d,%d:\t\t(%d,%d,%d)\t\t(%d,%d,%d)\n",i,j,force_index_i,x_row_index,y_row_index,z_row_index,x_column_index,y_column_index,z_column_index);
        // }

        a_matrix[x_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.x * Q1;
        a_matrix[x_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.x * Q2;
        a_matrix[x_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.x * Q3;
        a_matrix[y_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.y * Q1;
        a_matrix[y_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.y * Q2;
        a_matrix[y_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.y * Q3;
        a_matrix[z_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.z * Q1;
        a_matrix[z_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.z * Q2;
        a_matrix[z_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.z * Q3;

        if (force_index_i == 0)
        {
            //if (i == 0 && j == 1)
            //{
            //    printf("i=%d;j=%d\nBx=%f;By=%f;Bz=%f;D1=%f;QF=%f\n", i, j, b_vector[x_row_index], b_vector[y_row_index], b_vector[z_row_index], D1, QF);
            //}
            QF = TF1 * orientation_i.x + TF2 * orientation_i.y + TF3 * orientation_i.z;

            atomicAdd(&(b_vector[x_row_index]), -(D1 * orientation_i.x * QF));
            atomicAdd(&(b_vector[y_row_index]), -(D1 * orientation_i.y * QF));
            atomicAdd(&(b_vector[z_row_index]), -(D1 * orientation_i.z * QF));
            // b_vector[x_row_index] -= D1 * orientation_i.x * QF;
            // b_vector[y_row_index] -= D1 * orientation_i.y * QF;
            // b_vector[z_row_index] -= D1 * orientation_i.z * QF;
            //if (i == 0 && j == 1)
            //{
            //    printf("i=%d;j=%d\nBx=%f;By=%f;Bz=%f;D1=%f;QF=%f\n", i, j, b_vector[x_row_index], b_vector[y_row_index], b_vector[z_row_index], D1, QF);
            //}
        }

        for (force_index_j = 1; force_index_j < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index_j)
        {
            const fiberfloat gamma = 0.5f * (2.0f * (force_index_j + 1) + 1.0f) / (d + e - cc * lambda[force_index_j]);

            T11 = 0.0f;
            T22 = 0.0f;
            T33 = 0.0f;
            T12 = 0.0f;
            T13 = 0.0f;
            T23 = 0.0f;

            TF1 = 0.0f;
            TF2 = 0.0f;
            TF3 = 0.0f;

            for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
            {
                const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_i];
                const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_i + force_index_j * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
                T11 += quadrature_weight * G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T22 += quadrature_weight * G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T33 += quadrature_weight * G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T12 += quadrature_weight * G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T13 += quadrature_weight * G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T23 += quadrature_weight * G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;

                // if (i == 89 && j == 21 && force_index_i == 4 && force_index_j == 3 && quadrature_index_i==12)
                // {
                //     printf("%d,%d,%d,%d,%d,%f\n",i,j,force_index_i,force_index_j,quadrature_index_i,G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS]);
                // }

                if (force_index_i == 0)
                {
                    TF1 += quadrature_weight * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                    TF2 += quadrature_weight * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                    TF3 += quadrature_weight * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                }
            }

            Q1 = T11 * orientation_i.x + T12 * orientation_i.y + T13 * orientation_i.z;
            Q2 = T12 * orientation_i.x + T22 * orientation_i.y + T23 * orientation_i.z;
            Q3 = T13 * orientation_i.x + T23 * orientation_i.y + T33 * orientation_i.z;

            // if (i == 89 && j == 21 && force_index_i == 4 && force_index_j == 3)
            // {
            //     printf("i=%d;j=%d;force_index_i=%d;force_index_j=%d\nek=%f;gamma=%f;lambda=%f\n", i, j, force_index_i, force_index_j, eigen[force_index_j], gamma, lambda[force_index_j]);
            //     printf("i=%d;j=%d;force_index_i=%d;force_index_j=%d\nT11=%f;T22=%f;T33=%f;T12=%f;T13=%f;T23=%f\nQ1=%f\n", i, j, force_index_i, force_index_j, T11, T22, T33, T12, T13, T23, Q1);
            // }

            x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 0;
            y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 1;
            z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 2;

            x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 0;
            y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 1;
            z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 2;

            // if(x_row_index == 1344) {
            //     if(x_column_index == 327) {
            //         printf("xcol,%d,%d,%d,%d\n",i,j,force_index_i,force_index_j);
            //     }
            //     if(y_column_index == 327) {
            //         printf("ycol,%d,%d,%d,%d\n",i,j,force_index_i,force_index_j);
            //     }
            //     if(z_column_index == 327) {
            //         printf("zcol,%d,%d,%d,%d\n",i,j,force_index_i,force_index_j);
            //     }
            // }
            // if(y_row_index == 1344) {
            //     printf("yrow,%d,%d,%d,%d\n",i,j,force_index_i,force_index_j);
            // }
            // if(z_row_index == 1344) {
            //     printf("zrow,%d,%d,%d,%d\n",i,j,force_index_i,force_index_j);
            // }

            a_matrix[x_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T11 - eigen[force_index_j] * orientation_i.x * Q1);
            a_matrix[x_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T12 - eigen[force_index_j] * orientation_i.x * Q2);
            a_matrix[x_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T13 - eigen[force_index_j] * orientation_i.x * Q3);
            a_matrix[y_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T12 - eigen[force_index_j] * orientation_i.y * Q1);
            a_matrix[y_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T22 - eigen[force_index_j] * orientation_i.y * Q2);
            a_matrix[y_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T23 - eigen[force_index_j] * orientation_i.y * Q3);
            a_matrix[z_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T13 - eigen[force_index_j] * orientation_i.z * Q1);
            a_matrix[z_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T23 - eigen[force_index_j] * orientation_i.z * Q2);
            a_matrix[z_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T33 - eigen[force_index_j] * orientation_i.z * Q3);

            // if (i == 9 && j == 88 && force_index_i == 0 && force_index_j == 3)
            // {
            //     printf("%d,%d,%d,%d,%d,%d,%f\n", i, j, force_index_i, force_index_j, z_row_index, y_column_index, a_matrix[z_row_index + y_column_index * total_number_of_rows]);
            // }

            if (force_index_i == 0)
            {
                QF = TF1 * orientation_i.x + TF2 * orientation_i.y + TF3 * orientation_i.z;

                atomicAdd(&(b_vector[x_row_index]), -(gamma * (TF1 - eigen[force_index_j] * orientation_i.x * QF)));
                atomicAdd(&(b_vector[y_row_index]), -(gamma * (TF2 - eigen[force_index_j] * orientation_i.y * QF)));
                atomicAdd(&(b_vector[z_row_index]), -(gamma * (TF3 - eigen[force_index_j] * orientation_i.z * QF)));
            }
        }
    }
}

#endif // FIBERS_ASSEMBLE_SYSTEM_KERNEL_
