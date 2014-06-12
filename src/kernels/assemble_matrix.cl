kernel void assemble_matrix(const fiberuint number_of_fibers,
                            const fiberuint number_of_terms_in_force_expansion,
                            const global fiberfloat4 *positions,
                            const global fiberfloat4 *orientations,
                            global fiberfloat *a_matrix,
                            global fiberfloat *quadrature_points,
                            global fiberfloat *quadrature_weights,
                            global fiberfloat *legendre_polynomials) 
{
    size_t i = get_global_id(0);

    if (i >= number_of_fibers) return;

    const fiberfloat4 position_i = positions[i];
    const fiberfloat4 orientation_i = orientations[i];

    const fiberuint dimensions = 3;
    const fiberuint total_number_of_rows = number_of_fibers * number_of_terms_in_force_expansion * dimensions;

    for(fiberuint j = 0; j < number_of_fibers; ++j) 
    {
        const fiberfloat4 position_j = positions[j];
        const fiberfloat4 orientation_j = orientations[j];

        for(fiberuint force_index = 0; force_index < number_of_terms_in_force_expansion; ++force_index) 
        {
            fiberuint x_row_index = i * number_of_terms_in_force_expansion * dimensions + force_index * dimensions + 0;
            fiberuint y_row_index = i * number_of_terms_in_force_expansion * dimensions + force_index * dimensions + 1;
            fiberuint z_row_index = i * number_of_terms_in_force_expansion * dimensions + force_index * dimensions + 2;

            fiberuint x_column_index = j * number_of_terms_in_force_expansion * dimensions + force_index * dimensions + 0;
            fiberuint y_column_index = j * number_of_terms_in_force_expansion * dimensions + force_index * dimensions + 1;
            fiberuint z_column_index = j * number_of_terms_in_force_expansion * dimensions + force_index * dimensions + 2;

            // if(i==0) {
            //     printf("%d,%d,%d:\t\t(%d,%d,%d)\t\t(%d,%d,%d)\n",i,j,force_index,x_row_index,y_row_index,z_row_index,x_column_index,y_column_index,z_column_index);
            // }

            if(i == j) 
            {
                a_matrix[x_row_index + x_column_index * total_number_of_rows] = 1;
                a_matrix[y_row_index + y_column_index * total_number_of_rows] = 1;
                a_matrix[z_row_index + z_column_index * total_number_of_rows] = 1;
            }
            else 
            {
                a_matrix[x_row_index + x_column_index * total_number_of_rows] = 8;
                a_matrix[y_row_index + y_column_index * total_number_of_rows] = 8;
                a_matrix[z_row_index + z_column_index * total_number_of_rows] = 8;
            }
        }        
    }
}
