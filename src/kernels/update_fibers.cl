kernel void update_fibers(const global fiberfloat4 *previous_positions,
                                    const global fiberfloat4 *current_positions,
                                    global fiberfloat4 *next_positions,
                                    const global fiberfloat4 *previous_orientations,
                                    const global fiberfloat4 *current_orientations,
                                    global fiberfloat4 *next_orientations,
                                    const global fiberfloat4 *previous_translational_velocities,
                                    const global fiberfloat4 *current_translational_velocities,
                                    const global fiberfloat4 *previous_rotational_velocities,
                                    const global fiberfloat4 *current_rotational_velocities)
{
    size_t i = get_global_id(0);

    if (i >= NUMBER_OF_FIBERS) return;

    next_positions[i] = 4.0f / 3.0f * current_positions[i]
                        - 1.0f / 3.0f * previous_positions[i]
                        + 2.0f / 3.0f * TIMESTEP * (2.0f * current_translational_velocities[i] - previous_translational_velocities[i]);
    next_orientations[i] = 4.0f / 3.0f * current_orientations[i]
                        - 1.0f / 3.0f * previous_orientations[i]
                        + 2.0f / 3.0f * TIMESTEP * (2.0f * current_rotational_velocities[i] - previous_rotational_velocities[i]);
}