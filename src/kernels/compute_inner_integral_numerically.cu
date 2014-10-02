#ifndef FIBERS_COMPUTE_INNER_INTEGRAL_NUMERICALLY_KERNEL_
#define FIBERS_COMPUTE_INNER_INTEGRAL_NUMERICALLY_KERNEL_

#include "constants.cu"

__device__
void compute_G_numeric(
     const fiberfloat4 position_i,
     const fiberfloat4 orientation_i,
     const fiberfloat4 position_j,
     const fiberfloat4 orientation_j,
     const fiberuint force_index,
     const fiberfloat4 external_force,
     fiberfloat *G,
     fiberfloat *GF,
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

        G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;

        GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;

        fiberfloat4 position_on_fiber_i;
        position_on_fiber_i.x = position_i.x + quadrature_points[quadrature_index_i] * orientation_i.x;
        position_on_fiber_i.y = position_i.y + quadrature_points[quadrature_index_i] * orientation_i.y;
        position_on_fiber_i.z = position_i.z + quadrature_points[quadrature_index_i] * orientation_i.z;

        #pragma unroll
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

            const fiberfloat d1 = difference.x * difference.x;
            const fiberfloat d2 = difference.y * difference.y;
            const fiberfloat d3 = difference.z * difference.z;

            const fiberfloat invDistance = rsqrtf(d1 + d2 + d3);
            const fiberfloat invDistance3 = invDistance * invDistance * invDistance;
            const fiberfloat invDistance5 = invDistance3 * invDistance * invDistance;

            // equation 10
            // Note:    The outer product of a vector with itself is always a symmetric matrix
            //          so to save computation we only compute the upper triangle.
            // TODO calculation can be optimized (i.e. not dividing by distance, simpifing etc.)
            const fiberfloat K11 = invDistance
                                   + invDistance3 * d1
                                   + 2.0f * SLENDERNESS * SLENDERNESS * (invDistance3
                                           - 3.0f * invDistance5 * d1);
            const fiberfloat K22 = invDistance
                                   + invDistance3 * d2
                                   + 2.0f * SLENDERNESS * SLENDERNESS * (invDistance3
                                           - 3.0f * invDistance5 * d2);
            const fiberfloat K33 = invDistance
                                   + invDistance3 * d3
                                   + 2.0f * SLENDERNESS * SLENDERNESS * (invDistance3
                                           - 3.0f * invDistance5 * d3);
            const fiberfloat K12 = invDistance3 * difference.x * difference.y
                                   + 2.0f * SLENDERNESS * SLENDERNESS
                                   * -3.0f * invDistance5 * difference.x * difference.y;
            const fiberfloat K13 = invDistance3 * difference.x * difference.z
                                   + 2.0f * SLENDERNESS * SLENDERNESS
                                   * -3.0f * invDistance5 * difference.x * difference.z;
            const fiberfloat K23 = invDistance3 * difference.y * difference.z
                                   + 2.0f * SLENDERNESS * SLENDERNESS
                                   * -3.0f * invDistance5 * difference.y * difference.z;

            const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_j];
            const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_j + force_index * TOTAL_NUMBER_OF_QUADRATURE_POINTS];

            // @TEST http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html#1262
            // Kahan Summation Formula

            // if (quadrature_index_j == 0) {
            //     GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = quadrature_weight * (K11 * external_force.x + K12 * external_force.y + K13 * external_force.z);
            //     GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = quadrature_weight * (K12 * external_force.x + K22 * external_force.y + K23 * external_force.z);
            //     GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = quadrature_weight * (K13 * external_force.x + K23 * external_force.y + K33 * external_force.z);

            //     CGF0 = 0.0f;
            //     CGF1 = 0.0f;
            //     CGF2 = 0.0f;
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
