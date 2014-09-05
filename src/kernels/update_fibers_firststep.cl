kernel void update_fibers_firststep(const global fiberfloat4 *current_positions,
                                    global fiberfloat4 *next_positions,
                                    const global fiberfloat4 *current_orientations,
                                    global fiberfloat4 *next_orientations,
                                    const global fiberfloat4 *translational_velocities,
                                    const global fiberfloat4 *rotational_velocities)
{
    size_t i = get_global_id(0);

    if (i >= NUMBER_OF_FIBERS) return;

    next_positions[i] = current_positions[i] + TIMESTEP * translational_velocities[i];
    next_orientations[i] = normalize(current_orientations[i] + TIMESTEP * rotational_velocities[i]);
}