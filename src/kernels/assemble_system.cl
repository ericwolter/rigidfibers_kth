void *compute_G_analytic(fiberfloat4 position_i,
                         fiberfloat4 orientation_i,
                         fiberfloat4 position_j,
                         fiberfloat4 orientation_j,
                         fiberuint force_index,
                         fiberfloat4 external_force,
                         global fiberfloat *quadrature_points,
                         global fiberfloat *quadrature_weights,
                         global fiberfloat *legendre_polynomials,
                         fiberfloat *G,
                         fiberfloat *GF,
                         bool debug)
{
    const fiberuint k = force_index + 1;

    for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
    {
        const fiberfloat4 position_on_fiber_i = position_i + quadrature_points[quadrature_index_i] * orientation_i;

        // the difference vector between the current point on the fiber and the center
        // of the other fiber
        const fiberfloat4 R0 = position_on_fiber_i - position_j; // R_0
        const fiberfloat b = -2.0
                             * (R0.x * orientation_j.x
                                + R0.y * orientation_j.y
                                + R0.z * orientation_j.z);
        const fiberfloat c = R0.x * R0.x
                             + R0.y * R0.y
                             + R0.z * R0.z;

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

        fiberfloat I1[NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3];
        fiberfloat I3[NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3];
        fiberfloat I5[NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3];


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

        I1[0] = log(fabs(2.0 * s_upper + b + 2.0 * u_upper)) - log(fabs(2.0 * s_lower + b + 2.0 * u_lower));
        I1[1] = u_upper - u_lower + (-b / 2.0) * I1[0];

        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     printf("%d  %.16f\n",0,I1[0]);
        //     printf("%d  %.16f \n",1,I1[1]);
        // }

        I3[0] = (d < 1e-7) ?
                (-2.0 / pown(2.0 * s_upper + b, 2)) - (-2.0 / pown(2.0 * s_lower + b, 2)) :
                ((2.0 * s_upper + b) / (2.0 * d * u_upper)) - ((2.0 * s_lower + b) / (2.0 * d * u_lower));
        I3[1] = (-1.0 / u_upper) - (-1.0 / u_lower) - b / 2.0 * I3[0];

        I5[0] = (d < 1e-7) ?
                (-4.0 / pown(2.0 * s_upper + b, 4)) - (-4 / pown(2.0 * s_lower + b, 4)) :
                ((2.0 * s_upper + b) / (6.0 * d * pown(u_upper, 3))) - ((2.0 * s_lower + b) / (6.0 * d * pown(u_lower, 3))) + (2.0 / (3.0 * d)) * I3[0];
        I5[1] = (d < 1e-7) ?
                (-8.0 / (3.0 * pown(2.0 * s_upper + b, 3))) - (-8.0 / (3.0 * pown(2.0 * s_lower + b, 3))) - (b / 2.0) * I5[0] :
                (-(b * s_upper + 2.0 * c) / (6.0 * d * pown(u_upper, 3.0))) - (-(b * s_lower + 2.0 * c) / (6.0 * d * pown(u_lower, 3))) - (b / (3.0 * d)) * I3[0];

        if (c < climit)
        {
            for (fiberint n = 2; n < k + 3; ++n)
            {
                I1[n] = (pown(s_upper, n - 1) * u_upper) / n - (pown(s_lower, n - 1) * u_lower) / n
                        + ((1.0 - 2.0 * n) * b) / (2.0 * n) * I1[n - 1] - ((n - 1) * c) / n * I1[n - 2];

                fiberfloat test = (-pown(s_upper, n - 1) * u_upper) / (1 - (n + 1)) + (pown(s_lower, n - 1) * u_lower) / (1 - (n + 1))
                                  - (((0.5 - ((n + 1) - 1)) * b) / (1 - (n + 1))) * I1[n - 1] + ((((n + 1) - 2) * c) / (1 - (n + 1))) * I1[n - 2];

                // if (debug && quadrature_index_i == 12 && k == 5)
                // {
                //     printf("%d  %.16f \n",n,I1[n]);
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
                           ((pown(s_upper, n + 1) * u_upper) / (n + 2.0)) - ((pown(s_lower, n + 1) * u_lower) / (n + 2.0))
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
}

void *compute_G(fiberfloat4 position_i,
                fiberfloat4 orientation_i,
                fiberfloat4 position_j,
                fiberfloat4 orientation_j,
                fiberuint force_index,
                fiberfloat4 external_force,
                global fiberfloat *quadrature_points,
                global fiberfloat *quadrature_weights,
                global fiberfloat *legendre_polynomials,
                fiberfloat *G,
                fiberfloat *GF) // @TODO better names
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

kernel void assemble_system(const global fiberfloat4 *positions,
                            const global fiberfloat4 *orientations,
                            global fiberfloat *a_matrix,
                            global fiberfloat *b_vector,
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

    // @TODO Constant external force
    const fiberfloat4 external_force = 0.5f * (fiberfloat4)(0.0, 0.0, -1.0, 0.0);

    const fiberuint total_number_of_rows = NUMBER_OF_FIBERS * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS;

    fiberuint x_row_index;
    fiberuint y_row_index;
    fiberuint z_row_index;

    fiberuint x_column_index;
    fiberuint y_column_index;
    fiberuint z_column_index;

    fiberfloat lambda[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
    fiberfloat eigen[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
    lambda[0] = 2.0;
    eigen[0] = ((d - e - cc * lambda[0]) / 2.0) / (d - cc * lambda[0]);

    b_vector[i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + 0] = 0.0;
    b_vector[i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + 1] = 0.0;
    b_vector[i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + 2] = 0.0;

    for (fiberuint force_index = 1; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
    {
        lambda[force_index] = lambda[force_index - 1] + 2.0 / (force_index + 1);
        eigen[force_index] = ((d - e - cc * lambda[force_index]) / 2.0) / (d - cc * lambda[force_index]);

        x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 0;
        y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 1;
        z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 2;

        b_vector[x_row_index] = 0.0;
        b_vector[y_row_index] = 0.0;
        b_vector[z_row_index] = 0.0;
    }

    // if( i==0 ) {
    //     printf("i=%d;lambda=(%f,%f,%f,%f,%f);eigen=(%f,%f,%f,%f,%f)\n",i,lambda[0],lambda[1],lambda[2],lambda[3],lambda[4],eigen[0],eigen[1],eigen[2],eigen[3],eigen[4]);
    // }

    for (fiberuint j = 0; j < NUMBER_OF_FIBERS; ++j)
    {
        if (i == j)
        {
            // we only need to write the diagonals because the rest of the matrix
            // stays 0 throughout the simulation and is only initialised once
            // @TODO: why even do that on the GPU at all?
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

            // theta in equation 23
            fiberfloat T11 = 0.0;
            fiberfloat T22 = 0.0;
            fiberfloat T33 = 0.0;
            fiberfloat T12 = 0.0;
            fiberfloat T13 = 0.0;
            fiberfloat T23 = 0.0;

            fiberfloat TF1 = 0.0;
            fiberfloat TF2 = 0.0;
            fiberfloat TF3 = 0.0;
            fiberfloat QF;

            // @TODO combine computing G with the first iteration to calulate Theta(T11,...) for kk=1
            fiberfloat G[TOTAL_NUMBER_OF_QUADRATURE_POINTS * 6];
            fiberfloat GF[TOTAL_NUMBER_OF_QUADRATURE_POINTS * 3];
#ifdef USE_ANALYTICAL_INTEGRATION
            compute_G_analytic(position_i, orientation_i, position_j, orientation_j, force_index_i, external_force, quadrature_points, quadrature_weights, legendre_polynomials, G, GF, i == 89 && j == 21);
#else
            compute_G(position_i, orientation_i, position_j, orientation_j, force_index_i, external_force, quadrature_points, quadrature_weights, legendre_polynomials, G, GF);
#endif //USE_ANALYTICAL_INTEGRATION

            for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
            {
                const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_i];
                const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
                T11 += quadrature_weight * G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T22 += quadrature_weight * G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T33 += quadrature_weight * G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T12 += quadrature_weight * G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T13 += quadrature_weight * G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                T23 += quadrature_weight * G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;

                if (force_index_i == 0)
                {
                    TF1 += quadrature_weight * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                    TF2 += quadrature_weight * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                    TF3 += quadrature_weight * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                }
            }

            Q1 = T11 * orientation_i.x + T12 * orientation_i.y + T13 * orientation_i.z;
            Q2 = T12 * orientation_i.x + T22 * orientation_i.y + T23 * orientation_i.z;
            Q3 = T13 * orientation_i.x + T23 * orientation_i.y + T33 * orientation_i.z;

            //if (i == 0 && j == 1 && force_index_i == 0)
            //{
            //    printf("i=%d;j=%d;force_index_i=%d\nTF1=%f;TF2=%f;TF3=%f;QF=%f\n", i, j, force_index_i, TF1, TF2, TF3, QF);
            //}

            x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 0;
            y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 1;
            z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 2;

            x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 0;
            y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 1;
            z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 2;

            // if(i==0) {
            //     printf("%d,%d,%d:\t\t(%d,%d,%d)\t\t(%d,%d,%d)\n",i,j,force_index_i,x_row_index,y_row_index,z_row_index,x_column_index,y_column_index,z_column_index);
            // }

            a_matrix[x_row_index + x_column_index * total_number_of_rows] = D1 * orientation_i.x * Q1;
            a_matrix[x_row_index + y_column_index * total_number_of_rows] = D1 * orientation_i.x * Q2;
            a_matrix[x_row_index + z_column_index * total_number_of_rows] = D1 * orientation_i.x * Q3;
            a_matrix[y_row_index + x_column_index * total_number_of_rows] = D1 * orientation_i.y * Q1;
            a_matrix[y_row_index + y_column_index * total_number_of_rows] = D1 * orientation_i.y * Q2;
            a_matrix[y_row_index + z_column_index * total_number_of_rows] = D1 * orientation_i.y * Q3;
            a_matrix[z_row_index + x_column_index * total_number_of_rows] = D1 * orientation_i.z * Q1;
            a_matrix[z_row_index + y_column_index * total_number_of_rows] = D1 * orientation_i.z * Q2;
            a_matrix[z_row_index + z_column_index * total_number_of_rows] = D1 * orientation_i.z * Q3;

            if (force_index_i == 0)
            {
                //if (i == 0 && j == 1)
                //{
                //    printf("i=%d;j=%d\nBx=%f;By=%f;Bz=%f;D1=%f;QF=%f\n", i, j, b_vector[x_row_index], b_vector[y_row_index], b_vector[z_row_index], D1, QF);
                //}
                QF = TF1 * orientation_i.x + TF2 * orientation_i.y + TF3 * orientation_i.z;

                b_vector[x_row_index] -= D1 * orientation_i.x * QF;
                b_vector[y_row_index] -= D1 * orientation_i.y * QF;
                b_vector[z_row_index] -= D1 * orientation_i.z * QF;
                //if (i == 0 && j == 1)
                //{
                //    printf("i=%d;j=%d\nBx=%f;By=%f;Bz=%f;D1=%f;QF=%f\n", i, j, b_vector[x_row_index], b_vector[y_row_index], b_vector[z_row_index], D1, QF);
                //}
            }

            for (force_index_j = 1; force_index_j < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index_j)
            {
                fiberfloat gamma = 0.5 * (2.0 * (force_index_j + 1) + 1.0) / (d + e - cc * lambda[force_index_j]);

                fiberfloat T11 = 0.0;
                fiberfloat T22 = 0.0;
                fiberfloat T33 = 0.0;
                fiberfloat T12 = 0.0;
                fiberfloat T13 = 0.0;
                fiberfloat T23 = 0.0;

                fiberfloat TF1 = 0.0;
                fiberfloat TF2 = 0.0;
                fiberfloat TF3 = 0.0;

                for (fiberuint quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
                {
                    const fiberfloat quadrature_weight = quadrature_weights[quadrature_index_i];
                    const fiberfloat legendre_polynomial = legendre_polynomials[quadrature_index_i + force_index_j * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
                    T11 += quadrature_weight * G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                    T22 += quadrature_weight * G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                    T33 += quadrature_weight * G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                    T12 += quadrature_weight * G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                    T13 += quadrature_weight * G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                    T23 += quadrature_weight * G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;

                    // if (i == 89 && j == 21 && force_index_i == 4 && force_index_j == 3 && quadrature_index_i==12)
                    // {
                    //     printf("%d,%d,%d,%d,%d,%f\n",i,j,force_index_i,force_index_j,quadrature_index_i,G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS]);
                    // }

                    if (force_index_i == 0)
                    {
                        TF1 += quadrature_weight * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                        TF2 += quadrature_weight * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                        TF3 += quadrature_weight * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
                    }
                }

                Q1 = T11 * orientation_i.x + T12 * orientation_i.y + T13 * orientation_i.z;
                Q2 = T12 * orientation_i.x + T22 * orientation_i.y + T23 * orientation_i.z;
                Q3 = T13 * orientation_i.x + T23 * orientation_i.y + T33 * orientation_i.z;

                // if (i == 89 && j == 21 && force_index_i == 4 && force_index_j == 3)
                // {
                //     printf("i=%d;j=%d;force_index_i=%d;force_index_j=%d\nek=%f;gamma=%f;lambda=%f\n", i, j, force_index_i, force_index_j, eigen[force_index_j], gamma, lambda[force_index_j]);
                //     printf("i=%d;j=%d;force_index_i=%d;force_index_j=%d\nT11=%f;T22=%f;T33=%f;T12=%f;T13=%f;T23=%f\nQ1=%f\n", i, j, force_index_i, force_index_j, T11, T22, T33, T12, T13, T23, Q1);
                // }

                x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 0;
                y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 1;
                z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 2;

                x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 0;
                y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 1;
                z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_i * DIMENSIONS + 2;

                // if(x_row_index == 1344) {
                //     if(x_column_index == 327) {
                //         printf("xcol,%d,%d,%d,%d\n",i,j,force_index_i,force_index_j);
                //     }
                //     if(y_column_index == 327) {
                //         printf("ycol,%d,%d,%d,%d\n",i,j,force_index_i,force_index_j);
                //     }
                //     if(z_column_index == 327) {
                //         printf("zcol,%d,%d,%d,%d\n",i,j,force_index_i,force_index_j);
                //     }
                // }
                // if(y_row_index == 1344) {
                //     printf("yrow,%d,%d,%d,%d\n",i,j,force_index_i,force_index_j);
                // }
                // if(z_row_index == 1344) {
                //     printf("zrow,%d,%d,%d,%d\n",i,j,force_index_i,force_index_j);
                // }

                a_matrix[x_row_index + x_column_index * total_number_of_rows] = gamma * (T11 - eigen[force_index_j] * orientation_i.x * Q1);
                a_matrix[x_row_index + y_column_index * total_number_of_rows] = gamma * (T12 - eigen[force_index_j] * orientation_i.x * Q2);
                a_matrix[x_row_index + z_column_index * total_number_of_rows] = gamma * (T13 - eigen[force_index_j] * orientation_i.x * Q3);
                a_matrix[y_row_index + x_column_index * total_number_of_rows] = gamma * (T12 - eigen[force_index_j] * orientation_i.y * Q1);
                a_matrix[y_row_index + y_column_index * total_number_of_rows] = gamma * (T22 - eigen[force_index_j] * orientation_i.y * Q2);
                a_matrix[y_row_index + z_column_index * total_number_of_rows] = gamma * (T23 - eigen[force_index_j] * orientation_i.y * Q3);
                a_matrix[z_row_index + x_column_index * total_number_of_rows] = gamma * (T13 - eigen[force_index_j] * orientation_i.z * Q1);
                a_matrix[z_row_index + y_column_index * total_number_of_rows] = gamma * (T23 - eigen[force_index_j] * orientation_i.z * Q2);
                a_matrix[z_row_index + z_column_index * total_number_of_rows] = gamma * (T33 - eigen[force_index_j] * orientation_i.z * Q3);

                // if (i == 9 && j == 88 && force_index_i == 0 && force_index_j == 3)
                // {
                //     printf("%d,%d,%d,%d,%d,%d,%f\n", i, j, force_index_i, force_index_j, z_row_index, y_column_index, a_matrix[z_row_index + y_column_index * total_number_of_rows]);
                // }

                if (force_index_i == 0)
                {
                    QF = TF1 * orientation_i.x + TF2 * orientation_i.y + TF3 * orientation_i.z;

                    b_vector[x_row_index] -= gamma * (TF1 - eigen[force_index_j] * orientation_i.x * QF);
                    b_vector[y_row_index] -= gamma * (TF2 - eigen[force_index_j] * orientation_i.y * QF);
                    b_vector[z_row_index] -= gamma * (TF3 - eigen[force_index_j] * orientation_i.z * QF);
                }
            }
        }
    }
}
