#ifndef FIBERS_COMPUTE_INNER_INTEGRAL_NUMERICALLY_KERNEL_
#define FIBERS_COMPUTE_INNER_INTEGRAL_NUMERICALLY_KERNEL_

__device__
void compute_G_numeric(
     const fiberfloat4 position_i,
     const fiberfloat4 orientation_i,
     const fiberfloat4 position_j,
     const fiberfloat4 orientation_j,
     const fiberuint force_index,
     const fiberfloat4 external_force,
     const fiberfloat *quadrature_points,
     const fiberfloat *quadrature_weights,
     const fiberfloat *legendre_polynomials,
     fiberfloat *G,
     fiberfloat *GF,
     const fiberfloat SLENDERNESS,
     const fiberuint NUMBER_OF_TERMS_IN_FORCE_EXPANSION,
     const fiberuint TOTAL_NUMBER_OF_QUADRATURE_POINTS,
     const bool debug)
{
    for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
    {
        // fiberfloat CGF0;
        // fiberfloat CGF1;
        // fiberfloat CGF2;
        // fiberfloat YGF0;
        // fiberfloat YGF1;
        // fiberfloat YGF2;
        // fiberfloat TGF0;
        // fiberfloat TGF1;
        // fiberfloat TGF2;

        G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
        G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
        G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
        G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
        G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
        G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;

        GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
        GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
        GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;

        fiberfloat4 position_on_fiber_i;
        position_on_fiber_i.x = position_i.x + quadrature_points[quadrature_index_i] * orientation_i.x;
        position_on_fiber_i.y = position_i.y + quadrature_points[quadrature_index_i] * orientation_i.y;
        position_on_fiber_i.z = position_i.z + quadrature_points[quadrature_index_i] * orientation_i.z;

        for (fiberuint quadrature_index_j = 0; quadrature_index_j < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_j)
        {
            const fiberfloat quadrature_point = quadrature_points[quadrature_index_j];
            fiberfloat4 position_on_fiber_j;
            position_on_fiber_j.x = position_j.x + quadrature_point * orientation_j.x;
            position_on_fiber_j.y = position_j.y + quadrature_point * orientation_j.y;
            position_on_fiber_j.z = position_j.z + quadrature_point * orientation_j.z;

            fiberfloat4 difference;
            difference.x = position_on_fiber_i.x - position_on_fiber_j.x;
            difference.y = position_on_fiber_i.y - position_on_fiber_j.y;
            difference.z = position_on_fiber_i.z - position_on_fiber_j.z;

            const fiberfloat distance = sqrt(difference.x * difference.x + difference.y * difference.y + difference.z * difference.z);

            // equation 10
            // Note:    The outer product of a vector with itself is always a symmetric matrix
            //          so to save computation we only compute the upper triangle.
            // TODO calculation can be optimized (i.e. not dividing by distance, simpifing etc.)
            const fiberfloat K11 = 1.0 / distance
                                   + (1.0 / distance) * (difference.x / distance) * (difference.x / distance)
                                   + 2.0 * SLENDERNESS * SLENDERNESS * ((1.0 / (distance * distance * distance))
                                           - (3.0 / (distance * distance * distance)) * ((difference.x / distance) * (difference.x / distance)));
            const fiberfloat K22 = 1.0 / distance
                                   + (1.0 / distance) * (difference.y / distance) * (difference.y / distance)
                                   + 2.0 * SLENDERNESS * SLENDERNESS * ((1.0 / (distance * distance * distance))
                                           - (3.0 / (distance * distance * distance)) * ((difference.y / distance) * (difference.y / distance)));
            const fiberfloat K33 = 1.0 / distance
                                   + (1.0 / distance) * (difference.z / distance) * (difference.z / distance)
                                   + 2.0 * SLENDERNESS * SLENDERNESS * ((1.0 / (distance * distance * distance))
                                           - (3.0 / (distance * distance * distance)) * ((difference.z / distance) * (difference.z / distance)));
            const fiberfloat K12 = (1.0 / distance) * (difference.x / distance) * (difference.y / distance)
                                   + 2.0 * SLENDERNESS * SLENDERNESS
                                   * (-3.0 / (distance * distance * distance)) * (difference.x / distance) * (difference.y / distance);

            const fiberfloat K13 = (1.0 / distance) * (difference.x / distance) * (difference.z / distance)
                                   + 2.0 * SLENDERNESS * SLENDERNESS
                                   * (-3.0 / (distance * distance * distance)) * (difference.x / distance) * (difference.z / distance);

            const fiberfloat K23 = (1.0 / distance) * (difference.y / distance) * (difference.z / distance)
                                   + 2.0 * SLENDERNESS * SLENDERNESS
                                   * (-3.0 / (distance * distance * distance)) * (difference.y / distance) * (difference.z / distance);

            const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_j];
            const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_j + force_index * TOTAL_NUMBER_OF_QUADRATURE_POINTS];

            // @TEST http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html#1262
            // Kahan Summation Formula

            // if (quadrature_index_j == 0) {
            //     GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = quadrature_weight * (K11 * external_force.x + K12 * external_force.y + K13 * external_force.z);
            //     GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = quadrature_weight * (K12 * external_force.x + K22 * external_force.y + K23 * external_force.z);
            //     GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = quadrature_weight * (K13 * external_force.x + K23 * external_force.y + K33 * external_force.z);

            //     CGF0 = 0.0;
            //     CGF1 = 0.0;
            //     CGF2 = 0.0;
            // } else {
            //     YGF0 = quadrature_weight * (K11 * external_force.x + K12 * external_force.y + K13 * external_force.z) - CGF0;
            //     YGF1 = quadrature_weight * (K12 * external_force.x + K22 * external_force.y + K23 * external_force.z) - CGF1;
            //     YGF2 = quadrature_weight * (K13 * external_force.x + K23 * external_force.y + K33 * external_force.z) - CGF2;

            //     TGF0 = GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] + YGF0;
            //     TGF1 = GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] + YGF1;
            //     TGF2 = GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] + YGF2;

            //     CGF0 = (TGF0 - GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS]) -  YGF0;
            //     CGF1 = (TGF1 - GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS]) -  YGF1;
            //     CGF2 = (TGF2 - GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS]) -  YGF2;

            //     GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = TGF0;
            //     GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = TGF1;
            //     GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = TGF2;
            // }

            G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K11 * legendre_polynomial;
            G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K22 * legendre_polynomial;
            G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K33 * legendre_polynomial;
            G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K12 * legendre_polynomial;
            G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K13 * legendre_polynomial;
            G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K23 * legendre_polynomial;

            GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K11 * external_force.x + K12 * external_force.y + K13 * external_force.z);
            GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K12 * external_force.x + K22 * external_force.y + K23 * external_force.z);
            GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K13 * external_force.x + K23 * external_force.y + K33 * external_force.z);
        }
    }
}

#endif //FIBERS_COMPUTE_INNER_INTEGRAL_NUMERICALLY_KERNEL_
