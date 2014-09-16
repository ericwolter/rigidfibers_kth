#ifndef FIBERS_ASSEMBLE_SYSTEM_KERNEL_
#define FIBERS_ASSEMBLE_SYSTEM_KERNEL_

#include "../common.h"
#include "compute_inner_integral_analytically.cu"
#include "compute_inner_integral_numerically.cu"

#define DIMENSIONS 3

__global__ void assemble_system(
    const fiberfloat4 *positions,
    const fiberfloat4 *orientations,
    fiberfloat *a_matrix,
    fiberfloat *b_vector,
    const fiberfloat *quadrature_points,
    const fiberfloat *quadrature_weights,
    const fiberfloat *legendre_polynomials,
    const fiberuint NUMBER_OF_FIBERS,
    const fiberfloat SLENDERNESS,
    const fiberuint NUMBER_OF_TERMS_IN_FORCE_EXPANSION,
    const fiberuint TOTAL_NUMBER_OF_QUADRATURE_POINTS,
    const fiberint USE_ANALYTICAL_INTEGRATION
    )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= NUMBER_OF_FIBERS) return;

    const fiberfloat c  = log(SLENDERNESS * SLENDERNESS * M_E);
    const fiberfloat d  = -c;
    const fiberfloat e  = 2.0;
    const fiberfloat cc = 1.0;
    const fiberfloat D1 = 0.75 / (d - 2.0 * cc);

    const fiberfloat4 position_i = positions[i];
    const fiberfloat4 orientation_i = orientations[i];

    fiberfloat4 external_force;
    external_force.x = 0.5 * 0.0;
    external_force.y = 0.5 * 0.0;
    external_force.z = 0.5 * -1.0;

    const fiberuint total_number_of_rows = NUMBER_OF_FIBERS * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS;

    fiberuint x_row_index;
    fiberuint y_row_index;
    fiberuint z_row_index;

    fiberuint x_column_index;
    fiberuint y_column_index;
    fiberuint z_column_index;

    fiberfloat lambda[5];
    fiberfloat eigen[5];
    lambda[0] = 2.0;
    eigen[0] = ((d - e - cc * lambda[0]) / 2.0) / (d - cc * lambda[0]);

    b_vector[i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + 0] = 0.0;
    b_vector[i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + 1] = 0.0;
    b_vector[i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + 2] = 0.0;

    for (fiberuint force_index = 1; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
    {
        lambda[force_index] = lambda[force_index - 1] + 2.0 / (force_index + 1);
        eigen[force_index] = ((d - e - cc * lambda[force_index]) / 2.0) / (d - cc * lambda[force_index]);

        x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 0;
        y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 1;
        z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 2;

        b_vector[x_row_index] = 0.0;
        b_vector[y_row_index] = 0.0;
        b_vector[z_row_index] = 0.0;
    }

    for (fiberuint j = 0; j < NUMBER_OF_FIBERS; ++j)
    {
        if (i == j)
        {
            // we only need to write the diagonals because the rest of the matrix
            // stays 0 throughout the simulation and is only initialised once
            // @TODO: why even do that on the GPU at all?
            for (fiberuint force_index = 0; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
            {
                x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 0;
                y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 1;
                z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 2;

                x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 0;
                y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 1;
                z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 2;

                a_matrix[x_row_index + x_column_index * total_number_of_rows] = 1;
                a_matrix[y_row_index + y_column_index * total_number_of_rows] = 1;
                a_matrix[z_row_index + z_column_index * total_number_of_rows] = 1;
            }

            continue;
        }

        const fiberfloat4 position_j = positions[j];
        const fiberfloat4 orientation_j = orientations[j];

        for (fiberuint force_index_i = 0; force_index_i < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index_i)
        {
            fiberfloat Q1;
            fiberfloat Q2;
            fiberfloat Q3;

            fiberuint force_index_j = 0;

            // theta in equation 23
            fiberfloat T11 = 0.0;
            fiberfloat T22 = 0.0;
            fiberfloat T33 = 0.0;
            fiberfloat T12 = 0.0;
            fiberfloat T13 = 0.0;
            fiberfloat T23 = 0.0;

            fiberfloat TF1 = 0.0;
            fiberfloat TF2 = 0.0;
            fiberfloat TF3 = 0.0;
            fiberfloat QF;

            fiberfloat G[24 * 6];
            fiberfloat GF[24 * 3];

            if(USE_ANALYTICAL_INTEGRATION) {
                compute_G_analytic(position_i, orientation_i, position_j, orientation_j, force_index_i, external_force, quadrature_points, quadrature_weights, legendre_polynomials, G, GF, SLENDERNESS, NUMBER_OF_TERMS_IN_FORCE_EXPANSION, TOTAL_NUMBER_OF_QUADRATURE_POINTS, i == 89 && j == 21);
            } else {
                compute_G_numeric(position_i, orientation_i, position_j, orientation_j, force_index_i, external_force, quadrature_points, quadrature_weights, legendre_polynomials, G, GF, SLENDERNESS, NUMBER_OF_TERMS_IN_FORCE_EXPANSION, TOTAL_NUMBER_OF_QUADRATURE_POINTS, i == 89 && j == 21);
            }

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

            a_matrix[x_row_index + x_column_index * total_number_of_rows] = D1 * orientation_i.x * Q1;
            a_matrix[x_row_index + y_column_index * total_number_of_rows] = D1 * orientation_i.x * Q2;
            a_matrix[x_row_index + z_column_index * total_number_of_rows] = D1 * orientation_i.x * Q3;
            a_matrix[y_row_index + x_column_index * total_number_of_rows] = D1 * orientation_i.y * Q1;
            a_matrix[y_row_index + y_column_index * total_number_of_rows] = D1 * orientation_i.y * Q2;
            a_matrix[y_row_index + z_column_index * total_number_of_rows] = D1 * orientation_i.y * Q3;
            a_matrix[z_row_index + x_column_index * total_number_of_rows] = D1 * orientation_i.z * Q1;
            a_matrix[z_row_index + y_column_index * total_number_of_rows] = D1 * orientation_i.z * Q2;
            a_matrix[z_row_index + z_column_index * total_number_of_rows] = D1 * orientation_i.z * Q3;

            if (force_index_i == 0)
            {
                //if (i == 0 && j == 1)
                //{
                //    printf("i=%d;j=%d\nBx=%f;By=%f;Bz=%f;D1=%f;QF=%f\n", i, j, b_vector[x_row_index], b_vector[y_row_index], b_vector[z_row_index], D1, QF);
                //}
                QF = TF1 * orientation_i.x + TF2 * orientation_i.y + TF3 * orientation_i.z;

                b_vector[x_row_index] -= D1 * orientation_i.x * QF;
                b_vector[y_row_index] -= D1 * orientation_i.y * QF;
                b_vector[z_row_index] -= D1 * orientation_i.z * QF;
                //if (i == 0 && j == 1)
                //{
                //    printf("i=%d;j=%d\nBx=%f;By=%f;Bz=%f;D1=%f;QF=%f\n", i, j, b_vector[x_row_index], b_vector[y_row_index], b_vector[z_row_index], D1, QF);
                //}
            }

            for (force_index_j = 1; force_index_j < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index_j)
            {
                fiberfloat gamma = 0.5 * (2.0 * (force_index_j + 1) + 1.0) / (d + e - cc * lambda[force_index_j]);

                fiberfloat T11 = 0.0;
                fiberfloat T22 = 0.0;
                fiberfloat T33 = 0.0;
                fiberfloat T12 = 0.0;
                fiberfloat T13 = 0.0;
                fiberfloat T23 = 0.0;

                fiberfloat TF1 = 0.0;
                fiberfloat TF2 = 0.0;
                fiberfloat TF3 = 0.0;

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

                a_matrix[x_row_index + x_column_index * total_number_of_rows] = gamma * (T11 - eigen[force_index_j] * orientation_i.x * Q1);
                a_matrix[x_row_index + y_column_index * total_number_of_rows] = gamma * (T12 - eigen[force_index_j] * orientation_i.x * Q2);
                a_matrix[x_row_index + z_column_index * total_number_of_rows] = gamma * (T13 - eigen[force_index_j] * orientation_i.x * Q3);
                a_matrix[y_row_index + x_column_index * total_number_of_rows] = gamma * (T12 - eigen[force_index_j] * orientation_i.y * Q1);
                a_matrix[y_row_index + y_column_index * total_number_of_rows] = gamma * (T22 - eigen[force_index_j] * orientation_i.y * Q2);
                a_matrix[y_row_index + z_column_index * total_number_of_rows] = gamma * (T23 - eigen[force_index_j] * orientation_i.y * Q3);
                a_matrix[z_row_index + x_column_index * total_number_of_rows] = gamma * (T13 - eigen[force_index_j] * orientation_i.z * Q1);
                a_matrix[z_row_index + y_column_index * total_number_of_rows] = gamma * (T23 - eigen[force_index_j] * orientation_i.z * Q2);
                a_matrix[z_row_index + z_column_index * total_number_of_rows] = gamma * (T33 - eigen[force_index_j] * orientation_i.z * Q3);

                // if (i == 9 && j == 88 && force_index_i == 0 && force_index_j == 3)
                // {
                //     printf("%d,%d,%d,%d,%d,%d,%f\n", i, j, force_index_i, force_index_j, z_row_index, y_column_index, a_matrix[z_row_index + y_column_index * total_number_of_rows]);
                // }

                if (force_index_i == 0)
                {
                    QF = TF1 * orientation_i.x + TF2 * orientation_i.y + TF3 * orientation_i.z;

                    b_vector[x_row_index] -= gamma * (TF1 - eigen[force_index_j] * orientation_i.x * QF);
                    b_vector[y_row_index] -= gamma * (TF2 - eigen[force_index_j] * orientation_i.y * QF);
                    b_vector[z_row_index] -= gamma * (TF3 - eigen[force_index_j] * orientation_i.z * QF);
                }
            }
        }
    }

    // delete[] G;
    // delete[] GF;
}

#endif // FIBERS_ASSEMBLE_SYSTEM_KERNEL_
