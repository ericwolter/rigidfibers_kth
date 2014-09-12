#ifndef FIBERS_COMPUTE_INNER_INTEGRAL_ANALYTICALLY_KERNEL_
#define FIBERS_COMPUTE_INNER_INTEGRAL_ANALYTICALLY_KERNEL_

__device__
void compute_G_analytic(
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
    const fiberuint k = force_index + 1;

    fiberfloat *I1 = new fiberfloat[NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3];
    fiberfloat *I3 = new fiberfloat[NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3];
    fiberfloat *I5 = new fiberfloat[NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3];

    for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
    {
        fiberfloat4 position_on_fiber_i;
        position_on_fiber_i.x = position_i.x + quadrature_points[quadrature_index_i] * orientation_i.x;
        position_on_fiber_i.y = position_i.y + quadrature_points[quadrature_index_i] * orientation_i.y;
        position_on_fiber_i.z = position_i.z + quadrature_points[quadrature_index_i] * orientation_i.z;

        // the difference vector between the current point on the fiber and the center
        // of the other fiber
        fiberfloat4 R0;
        R0.x = position_on_fiber_i.x - position_j.x;
        R0.y = position_on_fiber_i.y - position_j.y;
        R0.z = position_on_fiber_i.z - position_j.z;
        const fiberfloat b = -2.0
                             * (R0.x * orientation_j.x
                                + R0.y * orientation_j.y
                                + R0.z * orientation_j.z);
        const fiberfloat c = R0.x * R0.x
                             + R0.y * R0.y
                             + R0.z * R0.z;
        //const fiberfloat c_test = ((R0.x * R0.y * R0.z) * (R0.x * R0.y * R0.z)) - 2 * R0.x * R0.y - 2 * R0.x * R0.z - 2 * R0.y * R0.z;

        // if fibers are too far apart we have numerical problems
        // so in order to minimize the effect we inverse the
        // the recursive direction
        // @TODO How/why does this help exactly?
        const fiberfloat climit = 10.0;

        const fiberfloat d = c - 0.25 * b * b;

        const fiberfloat s_upper = 1.0;
        const fiberfloat s_lower = -1.0;

        const fiberfloat u_upper = sqrt(s_upper * s_upper + b * s_upper + c);
        const fiberfloat u_lower = sqrt(s_lower * s_lower + b * s_lower + c);

        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     fiberfloat y = 0;
        //     const fiberfloat x = fabs(2.0 * s_upper + b + 2.0 * u_upper);

        //     for (int i = 0; i < 100; ++i)
        //     {
        //         y = y + 2 * ((x - exp(y))/(x + exp(y)));
        //         printf("%d, %f\n",i,y);
        //     }
        // }

        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     printf(" %.16f\n", 0, position_j.x);
        //     printf(" %.16f\n", 0, c);
        //     printf(" %.16f\n", 0, c_test);
        //     printf(" %.16f\n", 0, u_upper);
        //     printf(" %.16f\n", 0, 2.0 * s_upper + b + 2.0 * u_upper);
        //     printf(" %.16f\n", 0, fabs(2.0 * s_upper + b + 2.0 * u_upper));
        //     printf(" %.16f\n", 0, log(fabs(2.0 * s_upper + b + 2.0 * u_upper)));
        // }

        I1[0] = log(fabs(2.0 * s_upper + b + 2.0 * u_upper)) - log(fabs(2.0 * s_lower + b + 2.0 * u_lower));
        I1[1] = u_upper - u_lower + (-b / 2.0) * I1[0];

        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     printf("%d  %.16f\n", 0, I1[0]);
        //     printf("%d  %.16f \n", 1, I1[1]);
        // }

        I3[0] = (d < 1e-7) ?
                (-2.0 / pow(2.0 * s_upper + b, 2)) - (-2.0 / pow(2.0 * s_lower + b, 2)) :
                ((2.0 * s_upper + b) / (2.0 * d * u_upper)) - ((2.0 * s_lower + b) / (2.0 * d * u_lower));
        I3[1] = (-1.0 / u_upper) - (-1.0 / u_lower) - b / 2.0 * I3[0];

        I5[0] = (d < 1e-7) ?
                (-4.0 / pow(2.0 * s_upper + b, 4)) - (-4 / pow(2.0 * s_lower + b, 4)) :
                ((2.0 * s_upper + b) / (6.0 * d * pow(u_upper, 3))) - ((2.0 * s_lower + b) / (6.0 * d * pow(u_lower, 3))) + (2.0 / (3.0 * d)) * I3[0];
        I5[1] = (d < 1e-7) ?
                (-8.0 / (3.0 * pow(2.0 * s_upper + b, 3))) - (-8.0 / (3.0 * pow(2.0 * s_lower + b, 3))) - (b / 2.0) * I5[0] :
                (-(b * s_upper + 2.0 * c) / (6.0 * d * pow(u_upper, 3))) - (-(b * s_lower + 2.0 * c) / (6.0 * d * pow(u_lower, 3))) - (b / (3.0 * d)) * I3[0];

        if (c < climit)
        {
            for (fiberint n = 2; n < k + 3; ++n)
            {
                I1[n] = (pow(s_upper, n - 1) * u_upper) / n - (pow(s_lower, n - 1) * u_lower) / n
                        + ((1.0 - 2.0 * n) * b) / (2.0 * n) * I1[n - 1] - ((n - 1) * c) / n * I1[n - 2];

                fiberfloat test = (-pow(s_upper, n - 1) * u_upper) / (1 - (n + 1)) + (pow(s_lower, n - 1) * u_lower) / (1 - (n + 1))
                                  - (((0.5 - ((n + 1) - 1)) * b) / (1 - (n + 1))) * I1[n - 1] + ((((n + 1) - 2) * c) / (1 - (n + 1))) * I1[n - 2];

                // if (debug && quadrature_index_i == 12 && k == 5)
                // {
                //     printf("%d  %.16f \n", n, I1[n]);
                // }

                I3[n] = I1[n - 2] - b * I3[n - 1] - c * I3[n - 2];

                I5[n] = I3[n - 2] - b * I5[n - 1] - c * I5[n - 2];
            }
        }
        else
        {
            fiberfloat i1n2 = 0.0;
            fiberfloat i1n1 = 0.0;
            fiberfloat i1n0 = 0.0;

            fiberfloat i3n2 = 0.0;
            fiberfloat i3n1 = 0.0;
            fiberfloat i3n0 = 0.0;

            fiberfloat i5n2 = 0.0;
            fiberfloat i5n1 = 0.0;
            fiberfloat i5n0 = 0.0;

            // if (debug && quadrature_index_i == 0 && k == 5)
            // {
            //     printf("%d,%f\n", 29, i1n2);
            //     printf("%d,%f\n", 28, i1n1);
            // }

            for (fiberint n = 27; n > 1; --n)
            {
                i1n0 = (n + 2.0) / ((n + 1.0) * c) * (
                           ((pow(s_upper, n + 1) * u_upper) / (n + 2.0)) - ((pow(s_lower, n + 1) * u_lower) / (n + 2.0))
                           + ((1.0 - 2.0 * (n + 2.0)) / (2.0 * (n + 2.0))) * b * i1n1 - i1n2);
                i3n0 = 1.0 / c * (i1n0 - b * i3n1 - i3n2);
                i5n0 = 1.0 / c * (i3n0 - b * i5n1 - i5n2);

                // if (debug && quadrature_index_i == 0 && k == 5)
                // {
                //     printf("%d,%f,%f,%f,%f,%f\n", n, i1n0, i1n1, i1n2,
                //            ((pown(s_upper, n + 1) * u_upper) / (n + 2)) - ((pown(s_lower, n + 1) * u_lower) / (n + 2)),
                //            ((1.0 - 2.0 * (n + 2.0)) / (2.0 * (n + 2.0))) );
                // }

                if (n < NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3)
                {
                    I1[n] = i1n0;
                    I3[n] = i3n0;
                    I5[n] = i5n0;
                }

                i1n2 = i1n1;
                i1n1 = i1n0;

                i3n2 = i3n1;
                i3n1 = i3n0;

                i5n2 = i5n1;
                i5n1 = i5n0;
            }
        }

        fiberfloat L01;
        fiberfloat L03;
        fiberfloat L05;
        fiberfloat L13;
        fiberfloat L15;
        fiberfloat L23;
        fiberfloat L25;


        if (k == 0 || k == 1)
        {
            L01 = I1[k];
            L03 = I3[k];
            L05 = I5[k];
            L13 = I3[k + 1];
            L15 = I5[k + 1];
            L23 = I3[k + 2];
            L25 = I5[k + 2];
        }
        else if (k == 2)
        {
            L01 = 0.5 * (3.0 * I1[k] - I1[k - 2]);
            L03 = 0.5 * (3.0 * I3[k] - I3[k - 2]);
            L05 = 0.5 * (3.0 * I5[k] - I5[k - 2]);
            L13 = 0.5 * (3.0 * I3[k + 1] - I3[k - 1]);
            L15 = 0.5 * (3.0 * I5[k + 1] - I5[k - 1]);
            L23 = 0.5 * (3.0 * I3[k + 2] - I3[k]);
            L25 = 0.5 * (3.0 * I5[k + 2] - I5[k]);
        }
        else if (k == 3)
        {
            L01 = 0.5 * (5.0 * I1[k] - 3.0 * I1[k - 2]);
            L03 = 0.5 * (5.0 * I3[k] - 3.0 * I3[k - 2]);
            L05 = 0.5 * (5.0 * I5[k] - 3.0 * I5[k - 2]);
            L13 = 0.5 * (5.0 * I3[k + 1] - 3.0 * I3[k - 1]);
            L15 = 0.5 * (5.0 * I5[k + 1] - 3.0 * I5[k - 1]);
            L23 = 0.5 * (5.0 * I3[k + 2] - 3.0 * I3[k]);
            L25 = 0.5 * (5.0 * I5[k + 2] - 3.0 * I5[k]);
        }
        else if (k == 4)
        {
            L01 = 0.125 * (35.0 * I1[k] - 30.0 * I1[k - 2] + 3.0 * I1[k - 4]);
            L03 = 0.125 * (35.0 * I3[k] - 30.0 * I3[k - 2] + 3.0 * I3[k - 4]);
            L05 = 0.125 * (35.0 * I5[k] - 30.0 * I5[k - 2] + 3.0 * I5[k - 4]);
            L13 = 0.125 * (35.0 * I3[k + 1] - 30.0 * I3[k - 1] + 3.0 * I3[k - 3]);
            L15 = 0.125 * (35.0 * I5[k + 1] - 30.0 * I5[k - 1] + 3.0 * I5[k - 3]);
            L23 = 0.125 * (35.0 * I3[k + 2] - 30.0 * I3[k] + 3.0 * I3[k - 2]);
            L25 = 0.125 * (35.0 * I5[k + 2] - 30.0 * I5[k] + 3.0 * I5[k - 2]);
        }
        else if (k == 5)
        {
            L01 = 0.125 * (63.0 * I1[k] - 70.0 * I1[k - 2] + 15.0 * I1[k - 4]);
            L03 = 0.125 * (63.0 * I3[k] - 70.0 * I3[k - 2] + 15.0 * I3[k - 4]);
            L05 = 0.125 * (63.0 * I5[k] - 70.0 * I5[k - 2] + 15.0 * I5[k - 4]);
            L13 = 0.125 * (63.0 * I3[k + 1] - 70.0 * I3[k - 1] + 15.0 * I3[k - 3]);
            L15 = 0.125 * (63.0 * I5[k + 1] - 70.0 * I5[k - 1] + 15.0 * I5[k - 3]);
            L23 = 0.125 * (63.0 * I3[k + 2] - 70.0 * I3[k] + 15.0 * I3[k - 2]);
            L25 = 0.125 * (63.0 * I5[k + 2] - 70.0 * I5[k] + 15.0 * I5[k - 2]);
        }
        else if (k == 6)
        {
            L01 = 0.0625 * (231.0 * I1[k] - 315.0 * I1[k - 2] + 105.0 * I1[k - 4] - 5.0 * I1[k - 6]);
            L03 = 0.0625 * (231.0 * I3[k] - 315.0 * I3[k - 2] + 105.0 * I3[k - 4] - 5.0 * I3[k - 6]);
            L05 = 0.0625 * (231.0 * I5[k] - 315.0 * I5[k - 2] + 105.0 * I5[k - 4] - 5.0 * I5[k - 6]);
            L13 = 0.0625 * (231.0 * I3[k + 1] - 315.0 * I3[k - 1] + 105.0 * I3[k - 3] - 5.0 * I3[k - 5]);
            L15 = 0.0625 * (231.0 * I5[k + 1] - 315.0 * I5[k - 1] + 105.0 * I5[k - 3] - 5.0 * I5[k - 5]);
            L23 = 0.0625 * (231.0 * I3[k + 2] - 315.0 * I3[k] + 105.0 * I3[k - 2] - 5.0 * I3[k - 4]);
            L25 = 0.0625 * (231.0 * I5[k + 2] - 315.0 * I5[k] + 105.0 * I5[k - 2] - 5.0 * I5[k - 4]);
        }
        else if (k == 7)
        {
            L01 = 0.0625 * (429.0 * I1[k] - 693.0 * I1[k - 2] + 315.0 * I1[k - 4] - 35.0 * I1[k - 6]);
            L03 = 0.0625 * (429.0 * I3[k] - 693.0 * I3[k - 2] + 315.0 * I3[k - 4] - 35.0 * I3[k - 6]);
            L05 = 0.0625 * (429.0 * I5[k] - 693.0 * I5[k - 2] + 315.0 * I5[k - 4] - 35.0 * I5[k - 6]);
            L13 = 0.0625 * (429.0 * I3[k + 1] - 693.0 * I3[k - 1] + 315.0 * I3[k - 3] - 35.0 * I3[k - 5]);
            L15 = 0.0625 * (429.0 * I5[k + 1] - 693.0 * I5[k - 1] + 315.0 * I5[k - 3] - 35.0 * I5[k - 5]);
            L23 = 0.0625 * (429.0 * I3[k + 2] - 693.0 * I3[k] + 315.0 * I3[k - 2] - 35.0 * I3[k - 4]);
            L25 = 0.0625 * (429.0 * I5[k + 2] - 693.0 * I5[k] + 315.0 * I5[k - 2] - 35.0 * I5[k - 4]);
        }
        else if (k == 8)
        {
            L01 = 0.0078125 * (6435.0 * I1[k] - 12012.0 * I1[k - 2] + 6930.0 * I1[k - 4] - 1260.0 * I1[k - 6] + 35.0 * I1[k - 8]);
            L03 = 0.0078125 * (6435.0 * I3[k] - 12012.0 * I3[k - 2] + 6930.0 * I3[k - 4] - 1260.0 * I3[k - 6] + 35.0 * I3[k - 8]);
            L05 = 0.0078125 * (6435.0 * I5[k] - 12012.0 * I5[k - 2] + 6930.0 * I5[k - 4] - 1260.0 * I5[k - 6] + 35.0 * I5[k - 8]);
            L13 = 0.0078125 * (6435.0 * I3[k + 1] - 12012.0 * I3[k - 1] + 6930.0 * I3[k - 3] - 1260.0 * I3[k - 5] + 35.0 * I3[k - 7]);
            L15 = 0.0078125 * (6435.0 * I5[k + 1] - 12012.0 * I5[k - 1] + 6930.0 * I5[k - 3] - 1260.0 * I5[k - 5] + 35.0 * I5[k - 7]);
            L23 = 0.0078125 * (6435.0 * I3[k + 2] - 12012.0 * I3[k] + 6930.0 * I3[k - 2] - 1260.0 * I3[k - 4] + 35.0 * I3[k - 6]);
            L25 = 0.0078125 * (6435.0 * I5[k + 2] - 12012.0 * I5[k] + 6930.0 * I5[k - 2] - 1260.0 * I5[k - 4] + 35.0 * I5[k - 6]);
        }

        G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G11
            L01
            + L03 * (R0.x * R0.x)
            - L13 * (orientation_j.x * R0.x + R0.x * orientation_j.x)
            + L23 * (orientation_j.x * orientation_j.x)
            + 2.0 * SLENDERNESS * SLENDERNESS * (
                L03
                - 3.0 * L05 * (R0.x * R0.x)
                + 3.0 * L15 * (orientation_j.x * R0.x + R0.x * orientation_j.x)
                - 3.0 * L25 * (orientation_j.x * orientation_j.x)
            );
        G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G22
            L01
            + L03 * (R0.y * R0.y)
            - L13 * (orientation_j.y * R0.y + R0.y * orientation_j.y)
            + L23 * (orientation_j.y * orientation_j.y)
            + 2.0 * SLENDERNESS * SLENDERNESS * (
                L03
                - 3.0 * L05 * (R0.y * R0.y)
                + 3.0 * L15 * (orientation_j.y * R0.y + R0.y * orientation_j.y)
                - 3.0 * L25 * (orientation_j.y * orientation_j.y)
            );
        G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G33
            L01
            + L03 * (R0.z * R0.z)
            - L13 * (orientation_j.z * R0.z + R0.z * orientation_j.z)
            + L23 * (orientation_j.z * orientation_j.z)
            + 2.0 * SLENDERNESS * SLENDERNESS * (
                L03
                - 3.0 * L05 * (R0.z * R0.z)
                + 3.0 * L15 * (orientation_j.z * R0.z + R0.z * orientation_j.z)
                - 3.0 * L25 * (orientation_j.z * orientation_j.z)
            );
        G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G12
            L03 * (R0.x * R0.y)
            - L13 * (orientation_j.x * R0.y + R0.x * orientation_j.y)
            + L23 * (orientation_j.x * orientation_j.y)
            + 2.0 * SLENDERNESS * SLENDERNESS * (
                - 3.0 * L05 * (R0.x * R0.y)
                + 3.0 * L15 * (orientation_j.x * R0.y + R0.x * orientation_j.y)
                - 3.0 * L25 * (orientation_j.x * orientation_j.y)
            );
        G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G13
            L03 * (R0.x * R0.z)
            - L13 * (orientation_j.x * R0.z + R0.x * orientation_j.z)
            + L23 * (orientation_j.x * orientation_j.z)
            + 2.0 * SLENDERNESS * SLENDERNESS * (
                - 3.0 * L05 * (R0.x * R0.z)
                + 3.0 * L15 * (orientation_j.x * R0.z + R0.x * orientation_j.z)
                - 3.0 * L25 * (orientation_j.x * orientation_j.z)
            );
        G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G23
            L03 * (R0.y * R0.z)
            - L13 * (orientation_j.y * R0.z + R0.y * orientation_j.z)
            + L23 * (orientation_j.y * orientation_j.z)
            + 2.0 * SLENDERNESS * SLENDERNESS * (
                - 3.0 * L05 * (R0.y * R0.z)
                + 3.0 * L15 * (orientation_j.y * R0.z + R0.y * orientation_j.z)
                - 3.0 * L25 * (orientation_j.y * orientation_j.z)
            );

        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     printf("a:%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
        //             c,
        //             log(fabs(2.0 * s_upper + b + 2.0 * u_upper)),
        //            I1[0],
        //            I1[1],
        //            I1[2],
        //            I1[3],
        //            I1[4],
        //            I1[5],
        //            I1[6],
        //            I1[7]
        //           );
        // }
        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     printf("b:%f,%f,%f,%f,%f,%f,%f,%f\n",
        //         G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS],
        //            L01,
        //            L03,
        //            L05,
        //            L13,
        //            L15,
        //            L23,
        //            L25
        //           );
        // }

        // if (debug && quadrature_index_i==0 && force_index == 5) {
        //     printf("%f,%f,%f\n",
        //         orientation_j.x * R0.x + R0.x * orientation_j.x,
        //         L13,
        //         (orientation_j.x * R0.x + R0.x * orientation_j.x)*L13);
        // }

        if (force_index == 0)
        {
            L01 = I1[0];
            L03 = I3[0];
            L05 = I5[0];
            L13 = I3[0 + 1];
            L15 = I5[0 + 1];
            L23 = I3[0 + 2];
            L25 = I5[0 + 2];

            const fiberfloat G11 =
                L01
                + L03 * (R0.x * R0.x)
                - L13 * (orientation_j.x * R0.x + R0.x * orientation_j.x)
                + L23 * (orientation_j.x * orientation_j.x)
                + 2.0 * SLENDERNESS * SLENDERNESS * (
                    L03
                    - 3.0 * L05 * (R0.x * R0.x)
                    + 3.0 * L15 * (orientation_j.x * R0.x + R0.x * orientation_j.x)
                    - 3.0 * L25 * (orientation_j.x * orientation_j.x)
                );
            const fiberfloat G22 =
                L01
                + L03 * (R0.y * R0.y)
                - L13 * (orientation_j.y * R0.y + R0.y * orientation_j.y)
                + L23 * (orientation_j.y * orientation_j.y)
                + 2.0 * SLENDERNESS * SLENDERNESS * (
                    L03
                    - 3.0 * L05 * (R0.y * R0.y)
                    + 3.0 * L15 * (orientation_j.y * R0.y + R0.y * orientation_j.y)
                    - 3.0 * L25 * (orientation_j.y * orientation_j.y)
                );
            const fiberfloat G33 =
                L01
                + L03 * (R0.z * R0.z)
                - L13 * (orientation_j.z * R0.z + R0.z * orientation_j.z)
                + L23 * (orientation_j.z * orientation_j.z)
                + 2.0 * SLENDERNESS * SLENDERNESS * (
                    L03
                    - 3.0 * L05 * (R0.z * R0.z)
                    + 3.0 * L15 * (orientation_j.z * R0.z + R0.z * orientation_j.z)
                    - 3.0 * L25 * (orientation_j.z * orientation_j.z)
                );
            const fiberfloat G12 =
                L03 * (R0.x * R0.y)
                - L13 * (orientation_j.x * R0.y + R0.x * orientation_j.y)
                + L23 * (orientation_j.x * orientation_j.y)
                + 2.0 * SLENDERNESS * SLENDERNESS * (
                    - 3.0 * L05 * (R0.x * R0.y)
                    + 3.0 * L15 * (orientation_j.x * R0.y + R0.x * orientation_j.y)
                    - 3.0 * L25 * (orientation_j.x * orientation_j.y)
                );
            const fiberfloat G13 =
                L03 * (R0.x * R0.z)
                - L13 * (orientation_j.x * R0.z + R0.x * orientation_j.z)
                + L23 * (orientation_j.x * orientation_j.z)
                + 2.0 * SLENDERNESS * SLENDERNESS * (
                    - 3.0 * L05 * (R0.x * R0.z)
                    + 3.0 * L15 * (orientation_j.x * R0.z + R0.x * orientation_j.z)
                    - 3.0 * L25 * (orientation_j.x * orientation_j.z)
                );
            const fiberfloat G23 =
                L03 * (R0.y * R0.z)
                - L13 * (orientation_j.y * R0.z + R0.y * orientation_j.z)
                + L23 * (orientation_j.y * orientation_j.z)
                + 2.0 * SLENDERNESS * SLENDERNESS * (
                    - 3.0 * L05 * (R0.y * R0.z)
                    + 3.0 * L15 * (orientation_j.y * R0.z + R0.y * orientation_j.z)
                    - 3.0 * L25 * (orientation_j.y * orientation_j.z)
                );

            GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] =
                G11 * external_force.x + G12 * external_force.y + G13 * external_force.z;
            GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] =
                G12 * external_force.x + G22 * external_force.y + G23 * external_force.z;
            GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] =
                G13 * external_force.x + G23 * external_force.y + G33 * external_force.z;
        }
    }

    delete[] I1;
    delete[] I3;
    delete[] I5;
}

#endif //FIBERS_COMPUTE_INNER_INTEGRAL_ANALYTICALLY_KERNEL_
