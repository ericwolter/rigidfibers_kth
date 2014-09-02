const float8 *compute_G(fiberfloat4 position_i,
                            fiberfloat4 orientation_i,
                            fiberfloat4 position_j,
                            fiberfloat4 orientation_j,
                            fiberuint force_index,
                            global fiberfloat *quadrature_points,
                            global fiberfloat *quadrature_weights,
                            global fiberfloat *legendre_polynomials)
{
    float8 G[TOTAL_NUMBER_OF_QUADRATURE_POINTS];

    for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
    {
        const fiberfloat4 position_on_fiber_i = position_i + quadrature_points[quadrature_index_i] * orientation_i;

        float8 test;
        for (fiberuint quadrature_index_j = 0; quadrature_index_j < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_j)
        {
            const fiberfloat quadrature_point = quadrature_points[quadrature_index_j];
            const fiberfloat4 position_on_fiber_j = position_j + quadrature_point * orientation_j;
            const fiberfloat4 difference = position_on_fiber_i - position_on_fiber_j;

            const fiberfloat inv_distance = rsqrt(difference.x * difference.x + difference.y * difference.y + difference.z * difference.z);
            const fiberfloat inv_distance_pow_3 = inv_distance * inv_distance * inv_distance;
            const fiberfloat4 norm_difference = difference * inv_distance;

            const fiberfloat s = 2.0 * SLENDERNESS * SLENDERNESS * (inv_distance_pow_3);
            const fiberfloat r = 3.0 * inv_distance_pow_3;
            const fiberfloat sr = s * (-r);

            // equation 10
            // Note:    The outer product of a vector with itself is always a symmetric matrix
            //          so to save computation we only compute the upper triangle.
            // TODO calculation can be optimized (i.e. not dividing by distance, simpifing etc.)
            const fiberfloat K11 = inv_distance
                                   + (inv_distance) * (norm_difference.x * norm_difference.x)
                                   + s
                                           - (r) * (norm_difference.x * norm_difference.x);
            const fiberfloat K22 = inv_distance
                                   + (inv_distance) * (norm_difference.y * norm_difference.y)
                                   + s
                                           - (r) * (norm_difference.y * norm_difference.y);
            const fiberfloat K33 = inv_distance
                                   + (inv_distance) * (norm_difference.z * norm_difference.z)
                                   + s
                                           - (r) * (norm_difference.z * norm_difference.z);
            const fiberfloat K12 = norm_difference.x * norm_difference.y * ((inv_distance) + sr);
            const fiberfloat K13 = norm_difference.x * norm_difference.z * ((inv_distance) + sr);
            const fiberfloat K23 = norm_difference.y * norm_difference.z * ((inv_distance) + sr);

            const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_j];
            const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_j + force_index * TOTAL_NUMBER_OF_QUADRATURE_POINTS];

            const fiberfloat ql = quadrature_weight * legendre_polynomial;
            test += ql * (float8)(K11,K22,K33,K12,K13,K23,0,0);

            // G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K11 * legendre_polynomial;
            // G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K22 * legendre_polynomial;
            // G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K33 * legendre_polynomial;
            // G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K12 * legendre_polynomial;
            // G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K13 * legendre_polynomial;
            // G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K23 * legendre_polynomial;
        }

        G[quadrature_index_i] = test;
    }

    return G;
}

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

    fiberfloat lambda[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
    fiberfloat eigen[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
    lambda[0] = 2.0;
    eigen[0] = ((d - e - cc * lambda[0])/2.0) / (d - cc * lambda[0]);
    for (fiberuint force_index = 1; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
    {
        lambda[force_index] = lambda[force_index - 1] + 2.0 / (force_index + 1);
        eigen[force_index] = ((d - e - cc * lambda[force_index])/2.0) / (d - cc * lambda[force_index]);
    }

    // if( i==0 ) {
    //     printf("i=%d;lambda=(%f,%f,%f,%f,%f);eigen=(%f,%f,%f,%f,%f)\n",i,lambda[0],lambda[1],lambda[2],lambda[3],lambda[4],eigen[0],eigen[1],eigen[2],eigen[3],eigen[4]);
    // }

    for (fiberuint j = 0; j < NUMBER_OF_FIBERS; ++j)
    {
        fiberuint x_row_index;
        fiberuint y_row_index;
        fiberuint z_row_index;

        fiberuint x_column_index;
        fiberuint y_column_index;
        fiberuint z_column_index;

        if(i == j) 
        {
            // we only need to write the diagonals because the rest of the matrix
            // stays 0 throughout the simulation and is only initialised once
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

            float8 T;

            // TODO combine computing G with the first iteration to calulate Theta(T11,...) for kk=1
            float8 *G = compute_G(position_i, orientation_i, position_j, orientation_j, force_index_i, quadrature_points, quadrature_weights, legendre_polynomials);

            for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
            {
                const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_i];
                const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
                const fiberfloat ql = quadrature_weight * legendre_polynomial;

                T += ql * G[quadrature_index_i];
            }

            Q1 = T.s0 * orientation_i.x + T.s3 * orientation_i.y + T.s4 * orientation_i.z;
            Q2 = T.s3 * orientation_i.x + T.s1 * orientation_i.y + T.s5 * orientation_i.z;
            Q3 = T.s4 * orientation_i.x + T.s5 * orientation_i.y + T.s2 * orientation_i.z;

            // if (i == 0 && j == 1)
            // {
            //     printf("i=%d;j=%d;force_index_i=%d\nT11=%f;T22=%f;T33=%f;T12=%f;T13=%f;T23=%f\nQ1=%f\n", i, j, force_index_i, T11, T22, T33, T12, T13, T23, Q1);
            // }

            x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 0;
            y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 1;
            z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 2;

            x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 0;
            y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 1;
            z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 2;

            // if(i==0) {
            //     printf("%d,%d,%d:\t\t(%d,%d,%d)\t\t(%d,%d,%d)\n",i,j,force_index_i,x_row_index,y_row_index,z_row_index,x_column_index,y_column_index,z_column_index);
            // }

            const fiberfloat4 a = D1 * orientation_i * (fiberfloat4)(Q1, Q2, Q3, 0);

            a_matrix[x_row_index + x_column_index * total_number_of_rows] = a.x;
            a_matrix[x_row_index + y_column_index * total_number_of_rows] = a.y;
            a_matrix[x_row_index + z_column_index * total_number_of_rows] = a.z;
            a_matrix[y_row_index + x_column_index * total_number_of_rows] = a.x;
            a_matrix[y_row_index + y_column_index * total_number_of_rows] = a.y;
            a_matrix[y_row_index + z_column_index * total_number_of_rows] = a.z;
            a_matrix[z_row_index + x_column_index * total_number_of_rows] = a.x;
            a_matrix[z_row_index + y_column_index * total_number_of_rows] = a.y;
            a_matrix[z_row_index + z_column_index * total_number_of_rows] = a.z;

            for (force_index_j = 1; force_index_j < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index_j)
            {
                const fiberfloat gamma = 0.5 * (2.0 * (force_index_j + 1) + 1.0)/(d+e-cc*lambda[force_index_j]);

                T *= 0.0;

                for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
                {
                    const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_i];
                    const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_i + force_index_j * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
                    const fiberfloat ql = quadrature_weight * legendre_polynomial;

                    T += ql * G[quadrature_index_i];
                }

                Q1 = T.s0 * orientation_i.x + T.s3 * orientation_i.y + T.s4 * orientation_i.z;
                Q2 = T.s3 * orientation_i.x + T.s1 * orientation_i.y + T.s5 * orientation_i.z;
                Q3 = T.s4 * orientation_i.x + T.s5 * orientation_i.y + T.s2 * orientation_i.z;

                // if (i == 1 && j == 0 && force_index_i == 0 && force_index_j == 1)
                // {
                //     printf("i=%d;j=%d;force_index_i=%d;force_index_j=%d\nek=%f;gamma=%f;lambda=%f\n", i, j, force_index_i, force_index_j, eigen[force_index_j],gamma, lambda[force_index_j]);
                //     printf("i=%d;j=%d;force_index_i=%d;force_index_j=%d\nT11=%f;T22=%f;T33=%f;T12=%f;T13=%f;T23=%f\nQ1=%f\n", i, j, force_index_i, force_index_j, T11, T22, T33, T12, T13, T23, Q1);
                // }

                x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 0;
                y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 1;
                z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 2;

                x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 0;
                y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 1;
                z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 2;     

                const fiberfloat4 a = eigen[force_index_j] * orientation_i * (fiberfloat4)(Q1, Q2, Q3, 0);

                a_matrix[x_row_index + x_column_index * total_number_of_rows] = gamma * (T.s0 - a.x);
                a_matrix[x_row_index + y_column_index * total_number_of_rows] = gamma * (T.s3 - a.y);
                a_matrix[x_row_index + z_column_index * total_number_of_rows] = gamma * (T.s4 - a.z);
                a_matrix[y_row_index + x_column_index * total_number_of_rows] = gamma * (T.s3 - a.x);
                a_matrix[y_row_index + y_column_index * total_number_of_rows] = gamma * (T.s1 - a.y);
                a_matrix[y_row_index + z_column_index * total_number_of_rows] = gamma * (T.s5 - a.z);
                a_matrix[z_row_index + x_column_index * total_number_of_rows] = gamma * (T.s4 - a.x);
                a_matrix[z_row_index + y_column_index * total_number_of_rows] = gamma * (T.s5 - a.y);
                a_matrix[z_row_index + z_column_index * total_number_of_rows] = gamma * (T.s2 - a.z);
            }            
        }
    }
}
