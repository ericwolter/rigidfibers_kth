kernel void assemble_matrix(const global fiberfloat4 *positions,
                            const global fiberfloat4 *orientations,
                            global fiberfloat *a_matrix,
                            global fiberfloat *quadrature_points,
                            global fiberfloat *quadrature_weights,
                            global fiberfloat *legendre_polynomials)
{
    size_t i = get_global_id(0);

    if (i >= NUMBER_OF_FIBERS) return;

    const fiberfloat c  = log(SLENDERNESS * SLENDERNESS * M_E_F);
    const fiberfloat d  = -c;
    const fiberfloat e  = 2.0;
    const fiberfloat cc = 1.0;
    const fiberfloat D1 = 0.75 / (d - 2.0 * cc);

    const fiberfloat4 position_i = positions[i];
    const fiberfloat4 orientation_i = orientations[i];

    const fiberuint total_number_of_rows = NUMBER_OF_FIBERS * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS;

    for (fiberuint j = 0; j < NUMBER_OF_FIBERS; ++j)
    {
        const fiberfloat4 position_j = positions[j];
        const fiberfloat4 orientation_j = orientations[j];

        for (fiberuint force_index = 0; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
        {
            // theta in equation 23
            fiberfloat T11 = 0.0;
            fiberfloat T22 = 0.0;
            fiberfloat T33 = 0.0;
            fiberfloat T12 = 0.0;
            fiberfloat T13 = 0.0;
            fiberfloat T23 = 0.0;

            fiberfloat G[TOTAL_NUMBER_OF_QUADRATURE_POINTS * 6];

            for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
            {
                const fiberfloat4 position_on_fiber_i = position_i + quadrature_points[quadrature_index_i] * orientation_i;

                G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
                G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
                G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
                G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
                G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
                G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;

                for (fiberuint quadrature_index_j = 0; quadrature_index_j < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_j)
                {
                    const fiberfloat quadrature_point = quadrature_points[quadrature_index_j];
                    const fiberfloat4 position_on_fiber_j = position_j + quadrature_point * orientation_j;
                    const fiberfloat4 difference = position_on_fiber_i - position_on_fiber_j;
                    const fiberfloat distance = length(difference);

                    // equation 10
                    // Note:    The outer product of a vector with itself is always a symmetric matrix
                    //          so to save computation we only compute the upper triangle.
                    // TODO calculation can be optimized (i.e. not dividing by distance, simpifing etc.)
                    const fiberfloat K11 = 1.0 / distance
                                           + (1.0 / distance) * ((difference.x / distance) * (difference.x / distance))
                                           + 2.0 * SLENDERNESS * SLENDERNESS * ((1.0 / (distance * distance * distance))
                                                   - (3.0 / (distance * distance * distance)) * ((difference.x / distance) * (difference.x / distance)));
                    const fiberfloat K22 = 1.0 / distance
                                           + (1.0 / distance) * ((difference.y / distance) * (difference.y / distance))
                                           + 2.0 * SLENDERNESS * SLENDERNESS * ((1.0 / (distance * distance * distance))
                                                   - (3.0 / (distance * distance * distance)) * ((difference.y / distance) * (difference.y / distance)));
                    const fiberfloat K33 = 1.0 / distance
                                           + (1.0 / distance) * ((difference.z / distance) * (difference.z / distance))
                                           + 2.0 * SLENDERNESS * SLENDERNESS * ((1.0 / (distance * distance * distance))
                                                   - (3.0 / (distance * distance * distance)) * ((difference.y / distance) * (difference.y / distance)));
                    const fiberfloat K12 = (1.0 / distance) * (difference.x / distance) * (difference.y / distance) 
                                           + 2.0 * SLENDERNESS * SLENDERNESS * ((1.0 / (distance * distance * distance))
                                                   * (-3.0 / (distance * distance * distance)) * (difference.x / distance) * (difference.y / distance));
                    const fiberfloat K13 = (1.0 / distance) * (difference.x / distance) * (difference.z / distance) 
                                           + 2.0 * SLENDERNESS * SLENDERNESS * ((1.0 / (distance * distance * distance))
                                                   * (-3.0 / (distance * distance * distance)) * (difference.x / distance) * (difference.z / distance));
                    const fiberfloat K23 = (1.0 / distance) * (difference.y / distance) * (difference.z / distance) 
                                           + 2.0 * SLENDERNESS * SLENDERNESS * ((1.0 / (distance * distance * distance))
                                                   * (-3.0 / (distance * distance * distance)) * (difference.y / distance) * (difference.z / distance));

                    const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_j];
                    const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_j + force_index * TOTAL_NUMBER_OF_QUADRATURE_POINTS];

                    G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K11 * legendre_polynomial;
                    G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K22 * legendre_polynomial;
                    G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K33 * legendre_polynomial;
                    G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K12 * legendre_polynomial;
                    G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K13 * legendre_polynomial;
                    G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K23 * legendre_polynomial;
                }

                const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_i];
                const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
                T11 += quadrature_weight * G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T22 += quadrature_weight * G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T33 += quadrature_weight * G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T12 += quadrature_weight * G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T13 += quadrature_weight * G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T23 += quadrature_weight * G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            }

            const fiberfloat Q1 = T11 * orientation_i.x + T12 * orientation_i.y + T13 * orientation_i.z;
            const fiberfloat Q2 = T12 * orientation_i.x + T22 * orientation_i.y + T23 * orientation_i.z;
            const fiberfloat Q3 = T13 * orientation_i.x + T23 * orientation_i.y + T33 * orientation_i.z;

            if(i==0 && j==1) {
                printf("i=%d;j=%d;force_index=%d\nT11=%f;T22=%f;T33=%f;T12=%f;T13=%f;T23=%f\nQ1=%f\n",i,j,force_index,T11,T22,T33,T12,T13,T23,Q1);
            }

            fiberuint x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 0;
            fiberuint y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 1;
            fiberuint z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 2;

            fiberuint x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 0;
            fiberuint y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 1;
            fiberuint z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 2;

            // if(i==0) {
            //     printf("%d,%d,%d:\t\t(%d,%d,%d)\t\t(%d,%d,%d)\n",i,j,force_index,x_row_index,y_row_index,z_row_index,x_column_index,y_column_index,z_column_index);
            // }

            if (i == j)
            {
                a_matrix[x_row_index + x_column_index * total_number_of_rows] = 1;
                a_matrix[y_row_index + y_column_index * total_number_of_rows] = 1;
                a_matrix[z_row_index + z_column_index * total_number_of_rows] = 1;
            }
            else
            {
                a_matrix[x_row_index + x_column_index * total_number_of_rows] = D1 * orientation_i.x * Q1;
                a_matrix[x_row_index + y_column_index * total_number_of_rows] = D1 * orientation_i.x * Q2;
                a_matrix[x_row_index + z_column_index * total_number_of_rows] = D1 * orientation_i.x * Q3;
                a_matrix[y_row_index + x_column_index * total_number_of_rows] = D1 * orientation_i.y * Q1;
                a_matrix[y_row_index + y_column_index * total_number_of_rows] = D1 * orientation_i.y * Q2;
                a_matrix[y_row_index + z_column_index * total_number_of_rows] = D1 * orientation_i.y * Q3;
                a_matrix[z_row_index + x_column_index * total_number_of_rows] = D1 * orientation_i.z * Q1;
                a_matrix[z_row_index + y_column_index * total_number_of_rows] = D1 * orientation_i.z * Q2;
                a_matrix[z_row_index + z_column_index * total_number_of_rows] = D1 * orientation_i.z * Q3;
            }
        }
    }
}
