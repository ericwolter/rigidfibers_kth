void *compute_GV(fiberuint j,
                fiberfloat4 position_i,
                fiberfloat4 orientation_i,
                fiberfloat4 position_j,
                fiberfloat4 orientation_j,
                global fiberfloat *coefficients,
                fiberfloat4 external_force,
                global fiberfloat *quadrature_points,
                global fiberfloat *quadrature_weights,
                global fiberfloat *legendre_polynomials,
                fiberfloat *GF,
                bool debug) // @TODO better names
{
    for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
    {
        GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
        GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;
        GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0;

        const fiberfloat4 position_on_fiber_i = position_i + quadrature_points[quadrature_index_i] * orientation_i;

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

            fiberfloat4 force_on_fiber_j = 0.5f * external_force;
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

            if (debug && quadrature_index_i == 15 - 1)
            {
                printf("     %f\n", K11);
            }
        }
    }
}

kernel void update_velocities(const global fiberfloat4 *positions,
                             const global fiberfloat4 *orientations,
                             global fiberfloat *coefficients,
                             global fiberfloat4 *translational_velocities,
                             global fiberfloat4 *rotational_velocities,
                             global fiberfloat *quadrature_points,
                             global fiberfloat *quadrature_weights,
                             global fiberfloat *legendre_polynomials)
{
    size_t i = get_global_id(0);

    if (i >= NUMBER_OF_FIBERS) return;

    const fiberfloat c  = log(SLENDERNESS * SLENDERNESS * M_E_F);
    const fiberfloat d  = -c;
    const fiberfloat e  = 2.0;

    const fiberfloat4 position_i = positions[i];
    const fiberfloat4 orientation_i = orientations[i];

    // @TODO Constant external force
    const fiberfloat4 external_force = (fiberfloat4)(0.0, 0.0, -1.0, 0.0);

    const fiberfloat4 oriented_force = (fiberfloat4)(
                                           orientation_i.x * orientation_i.x * external_force.x + orientation_i.x * orientation_i.y * external_force.y + orientation_i.x * orientation_i.z * external_force.z,
                                           orientation_i.x * orientation_i.y * external_force.x + orientation_i.y * orientation_i.y * external_force.y + orientation_i.y * orientation_i.z * external_force.z,
                                           orientation_i.x * orientation_i.z * external_force.x + orientation_i.y * orientation_i.z * external_force.y + orientation_i.z * orientation_i.z * external_force.z,
                                           0.0
                                       );

    translational_velocities[i] = 0.5f * ((d + 2) * external_force + (d - 2) * oriented_force);
    rotational_velocities[i] = (fiberfloat4)(0.0);

    for (fiberuint j = 0; j < NUMBER_OF_FIBERS; ++j)
    {
        if (i == j) continue;

        const fiberfloat4 position_j = positions[j];
        const fiberfloat4 orientation_j = orientations[j];

        fiberfloat GF[TOTAL_NUMBER_OF_QUADRATURE_POINTS * 3];
        compute_GV(j, position_i, orientation_i, position_j, orientation_j, coefficients, external_force, quadrature_points, quadrature_weights, legendre_polynomials, GF, false);

        fiberfloat TF1A0 = 0.0;
        fiberfloat TF2A0 = 0.0;
        fiberfloat TF3A0 = 0.0;

        fiberfloat TF1A1 = 0.0;
        fiberfloat TF2A1 = 0.0;
        fiberfloat TF3A1 = 0.0;

        for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
        {
            const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_i];
            const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];

            TF1A0 += quadrature_weight * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
            TF2A0 += quadrature_weight * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
            TF3A0 += quadrature_weight * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];

            TF1A1 += quadrature_weight * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            TF2A1 += quadrature_weight * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            TF3A1 += quadrature_weight * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
        }

        translational_velocities[i] += 0.5f * (fiberfloat4)(TF1A0, TF2A0, TF3A0, 0.0);
        rotational_velocities[i] += 1.5f * (fiberfloat4)(
                                        TF1A1 - (orientation_i.x * orientation_i.x * TF1A1 + orientation_i.x * orientation_i.y * TF2A1 + orientation_i.x * orientation_i.z * TF3A1),
                                        TF2A1 - (orientation_i.x * orientation_i.y * TF1A1 + orientation_i.y * orientation_i.y * TF2A1 + orientation_i.y * orientation_i.z * TF3A1),
                                        TF3A1 - (orientation_i.x * orientation_i.z * TF1A1 + orientation_i.y * orientation_i.z * TF2A1 + orientation_i.z * orientation_i.z * TF3A1),
                                        0.0
                                    );
    }

    translational_velocities[i] /= d;
    rotational_velocities[i] /= d;
}