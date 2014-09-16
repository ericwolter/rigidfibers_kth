#ifndef FIBERS_UPDATE_VELOCITIES_KERNEL_
#define FIBERS_UPDATE_VELOCITIES_KERNEL_

#include "constants.cu"

__device__
    void compute_GV(const fiberuint j,
                const fiberfloat4 position_i,
                const fiberfloat4 orientation_i,
                const fiberfloat4 position_j,
                const fiberfloat4 orientation_j,
                const fiberfloat *coefficients,
                const fiberfloat4 external_force,
                fiberfloat *GF
                ) // @TODO better names
{
    for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
    {
        GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;

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

            const fiberfloat d1 = difference.x * difference.x;
            const fiberfloat d2 = difference.x * difference.x;
            const fiberfloat d3 = difference.x * difference.x;

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
                                           - 3.0f * invDistance5 * difference.x * difference.x);
            const fiberfloat K22 = invDistance
                                   + invDistance3 * d2
                                   + 2.0f * SLENDERNESS * SLENDERNESS * (invDistance3
                                           - 3.0f * invDistance5 * difference.y * difference.y);
            const fiberfloat K33 = invDistance
                                   + invDistance3 * d3
                                   + 2.0f * SLENDERNESS * SLENDERNESS * (invDistance3
                                           - 3.0f * invDistance5 * difference.z * difference.z);
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

            fiberfloat4 force_on_fiber_j;
            force_on_fiber_j.x = 0.5f * external_force.x;
            force_on_fiber_j.y = 0.5f * external_force.y;
            force_on_fiber_j.z = 0.5f * external_force.z;

            for (fiberuint force_index = 0; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
            {
                const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_j + force_index * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
                fiberuint x_row_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 0;
                fiberuint y_row_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 1;
                fiberuint z_row_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 2;

                force_on_fiber_j.x += coefficients[x_row_index] * legendre_polynomial;
                force_on_fiber_j.y += coefficients[y_row_index] * legendre_polynomial;
                force_on_fiber_j.z += coefficients[z_row_index] * legendre_polynomial;
            }

            GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K11 * force_on_fiber_j.x + K12 * force_on_fiber_j.y + K13 * force_on_fiber_j.z);
            GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K12 * force_on_fiber_j.x + K22 * force_on_fiber_j.y + K23 * force_on_fiber_j.z);
            GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K13 * force_on_fiber_j.x + K23 * force_on_fiber_j.y + K33 * force_on_fiber_j.z);
        }
    }
}

__global__ void update_velocities(
    const fiberfloat4 *positions,
    const fiberfloat4 *orientations,
    const fiberfloat *coefficients,
    fiberfloat4 *translational_velocities,
    fiberfloat4 *rotational_velocities
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= NUMBER_OF_FIBERS) return;

    const fiberfloat c  = logf(SLENDERNESS * SLENDERNESS * M_E);
    const fiberfloat d  = -c;

    const fiberfloat4 position_i = positions[i];
    const fiberfloat4 orientation_i = orientations[i];

    // @TODO Constant external force
    fiberfloat4 external_force;
    external_force.x = 0.0f;
    external_force.y = 0.0f;
    external_force.z = -1.0f;

    fiberfloat4 oriented_force;
    oriented_force.x = orientation_i.x * orientation_i.x * external_force.x + orientation_i.x * orientation_i.y * external_force.y + orientation_i.x * orientation_i.z * external_force.z;
    oriented_force.y = orientation_i.x * orientation_i.y * external_force.x + orientation_i.y * orientation_i.y * external_force.y + orientation_i.y * orientation_i.z * external_force.z;
    oriented_force.z = orientation_i.x * orientation_i.z * external_force.x + orientation_i.y * orientation_i.z * external_force.y + orientation_i.z * orientation_i.z * external_force.z;

    translational_velocities[i].x = 0.5f * ((d + 2.0f) * external_force.x + (d - 2.0f) * oriented_force.x);
    translational_velocities[i].y = 0.5f * ((d + 2.0f) * external_force.y + (d - 2.0f) * oriented_force.y);
    translational_velocities[i].z = 0.5f * ((d + 2.0f) * external_force.z + (d - 2.0f) * oriented_force.z);

    rotational_velocities[i].x = 0.0f;
    rotational_velocities[i].y = 0.0f;
    rotational_velocities[i].z = 0.0f;

    for (fiberuint j = 0; j < NUMBER_OF_FIBERS; ++j)
    {
        if (i == j) continue;

        const fiberfloat4 position_j = positions[j];
        const fiberfloat4 orientation_j = orientations[j];

        fiberfloat GF[TOTAL_NUMBER_OF_QUADRATURE_POINTS * 3];
        compute_GV(j, position_i, orientation_i, position_j, orientation_j, coefficients, external_force, GF);

        fiberfloat TF1A0 = 0.0f;
        fiberfloat TF2A0 = 0.0f;
        fiberfloat TF3A0 = 0.0f;

        fiberfloat TF1A1 = 0.0f;
        fiberfloat TF2A1 = 0.0f;
        fiberfloat TF3A1 = 0.0f;

        for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
        {
            const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_i];
            const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
            const fiberfloat weighted_polynomial = quadrature_weight * legendre_polynomial;

            TF1A0 += quadrature_weight * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
            TF2A0 += quadrature_weight * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
            TF3A0 += quadrature_weight * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];

            TF1A1 += weighted_polynomial * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
            TF2A1 += weighted_polynomial * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
            TF3A1 += weighted_polynomial * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
        }

        translational_velocities[i].x += TF1A0;
        translational_velocities[i].y += TF2A0;
        translational_velocities[i].z += TF3A0;

        const fiberfloat t1 = orientation_i.x * TF1A1;
        const fiberfloat t2 = orientation_i.y * TF2A1;
        const fiberfloat t3 = orientation_i.z * TF3A1;

        rotational_velocities[i].x += TF1A1 - (orientation_i.x * t1 + orientation_i.x * t2 + orientation_i.x * t3);
        rotational_velocities[i].y += TF2A1 - (orientation_i.y * t1 + orientation_i.y * t2 + orientation_i.y * t3);
        rotational_velocities[i].z += TF3A1 - (orientation_i.z * t1 + orientation_i.z * t2 + orientation_i.z * t3);
    }

    translational_velocities[i].x *= (0.5f / d);
    translational_velocities[i].y *= (0.5f / d);
    translational_velocities[i].z *= (0.5f / d);

    rotational_velocities[i].x *= (1.5f / d);
    rotational_velocities[i].y *= (1.5f / d);
    rotational_velocities[i].z *= (1.5f / d);
}

#endif //FIBERS_UPDATE_VELOCITIES_KERNEL_
