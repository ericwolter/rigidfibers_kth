#ifndef FIBERS_UPDATE_FIBERS_KERNEL_
#define FIBERS_UPDATE_FIBERS_KERNEL_

__global__
void update_fibers(
    const fiberfloat4 *previous_positions,
    const fiberfloat4 *current_positions,
    fiberfloat4 *next_positions,
    const fiberfloat4 *previous_orientations,
    const fiberfloat4 *current_orientations,
    fiberfloat4 *next_orientations,
    const fiberfloat4 *previous_translational_velocities,
    const fiberfloat4 *current_translational_velocities,
    const fiberfloat4 *previous_rotational_velocities,
    const fiberfloat4 *current_rotational_velocities,
    const fiberuint NUMBER_OF_FIBERS,
    const fiberfloat TIMESTEP
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= NUMBER_OF_FIBERS) return;

    next_positions[i].x = 4.0 / 3.0 * current_positions[i].x
                        - 1.0 / 3.0 * previous_positions[i].x
                        + 2.0 / 3.0 * TIMESTEP * (2.0 * current_translational_velocities[i].x - previous_translational_velocities[i].x);
    next_positions[i].y = 4.0 / 3.0 * current_positions[i].y
                        - 1.0 / 3.0 * previous_positions[i].y
                        + 2.0 / 3.0 * TIMESTEP * (2.0 * current_translational_velocities[i].y - previous_translational_velocities[i].y);
    next_positions[i].z = 4.0 / 3.0 * current_positions[i].z
                        - 1.0 / 3.0 * previous_positions[i].z
                        + 2.0 / 3.0 * TIMESTEP * (2.0 * current_translational_velocities[i].z - previous_translational_velocities[i].z);

    next_orientations[i].x = 4.0 / 3.0 * current_orientations[i].x
                        - 1.0 / 3.0 * previous_orientations[i].x
                        + 2.0 / 3.0 * TIMESTEP * (2.0 * current_rotational_velocities[i].x - previous_rotational_velocities[i].x);
    next_orientations[i].y = 4.0 / 3.0 * current_orientations[i].y
                        - 1.0 / 3.0 * previous_orientations[i].y
                        + 2.0 / 3.0 * TIMESTEP * (2.0 * current_rotational_velocities[i].y - previous_rotational_velocities[i].y);
    next_orientations[i].z = 4.0 / 3.0 * current_orientations[i].z
                        - 1.0 / 3.0 * previous_orientations[i].z
                        + 2.0 / 3.0 * TIMESTEP * (2.0 * current_rotational_velocities[i].z - previous_rotational_velocities[i].z);

    fiberfloat invLen = 1.0 / sqrt(next_orientations[i].x * next_orientations[i].x 
        + next_orientations[i].y * next_orientations[i].y 
        + next_orientations[i].z * next_orientations[i].z);

    next_orientations[i].x *= invLen;
    next_orientations[i].y *= invLen;
    next_orientations[i].z *= invLen;
}

#endif //FIBERS_UPDATE_FIBERS_KERNEL_