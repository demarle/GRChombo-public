/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// code included in ApparentHorizon.impl.hpp::diff to calculate stencils

{
    const bool periodic[2] = { m_periodic_u, m_periodic_v };

    int coord[2] = { u, v };
    const int N[2] = { m_num_global_u, m_num_global_v };

    int d_start[2];
    double *const d_weights[2] = { out.du_weights, out.dv_weights };

    int dd_start[2];
    double *const dd_weights[2] = { out.dudu_weights, out.dvdv_weights };

    for (int i = 0; i < 2; ++i) {

        // Non-periodic boundaries
        if (!periodic[i] && (coord[i] == 0 || coord[i] == 1 || coord[i] == N[i] - 2 || coord[i] == N[i] - 1)) {

            if (coord[i] == 0) {

                d_start[i] = 0;

                d_weights[i][0] = -2.083333333333333333333;
                d_weights[i][1] = +4.000000000000000000000;
                d_weights[i][2] = -3.000000000000000000000;
                d_weights[i][3] = +1.333333333333333333333;
                d_weights[i][4] = -0.250000000000000000000;

                dd_start[i] = 0;

                dd_weights[i][0] = +3.75000000000000000000;
                dd_weights[i][1] = -12.8333333333333333333;
                dd_weights[i][2] = +17.8333333333333333333;
                dd_weights[i][3] = -13.0000000000000000000;
                dd_weights[i][4] = +5.08333333333333333333;
                dd_weights[i][5] = -0.83333333333333333333;
                
            }
            else if (coord[i] == 1) {

                d_start[i] = -1;

                d_weights[i][0] = -0.250000000000000000000;
                d_weights[i][1] = -0.833333333333333333333;
                d_weights[i][2] = +1.500000000000000000000;
                d_weights[i][3] = -0.500000000000000000000;
                d_weights[i][4] = +0.083333333333333333333;

                dd_start[i] = -1;

                dd_weights[i][0] = +0.83333333333333333333;
                dd_weights[i][1] = -1.25000000000000000000;
                dd_weights[i][2] = -0.33333333333333333333;
                dd_weights[i][3] = +1.16666666666666666667;
                dd_weights[i][4] = -0.50000000000000000000;
                dd_weights[i][5] = +0.08333333333333333333;
                
            }
            else if (coord[i] == N[i] - 2) {

                d_start[i] = -3;

                d_weights[i][0] = -0.083333333333333333333;
                d_weights[i][1] = +0.500000000000000000000;
                d_weights[i][2] = -1.500000000000000000000;
                d_weights[i][3] = +0.833333333333333333333;
                d_weights[i][4] = +0.250000000000000000000;

                dd_start[i] = -4;

                dd_weights[i][0] = +0.08333333333333333333;
                dd_weights[i][1] = -0.50000000000000000000;
                dd_weights[i][2] = +1.16666666666666666667;
                dd_weights[i][3] = -0.33333333333333333333;
                dd_weights[i][4] = -1.25000000000000000000;
                dd_weights[i][5] = +0.83333333333333333333;
                
            }
            else if (coord[i] == N[i] - 1) {

                d_start[i] = -4;

                d_weights[i][0] = +0.250000000000000000000;
                d_weights[i][1] = -1.333333333333333333333;
                d_weights[i][2] = +3.000000000000000000000;
                d_weights[i][3] = -4.000000000000000000000;
                d_weights[i][4] = +2.083333333333333333333;

                dd_start[i] = -5;

                dd_weights[i][0] = -0.83333333333333333333;
                dd_weights[i][1] = +5.08333333333333333333;
                dd_weights[i][2] = -13.0000000000000000000;
                dd_weights[i][3] = +17.8333333333333333333;
                dd_weights[i][4] = -12.8333333333333333333;
                dd_weights[i][5] = +3.75000000000000000000;
                
            }
        }
        // Standard points
        else {

            d_start[i] = -2;

            d_weights[i][0] = +0.083333333333333333333;
            d_weights[i][1] = -0.666666666666666666666;
            d_weights[i][2] = 0;
            d_weights[i][3] = +0.666666666666666666666;
            d_weights[i][4] = -0.083333333333333333333;

            if (coord[i] == N[i] - 3)
            {
                dd_start[i] = -3;

                dd_weights[i][0] = 0;
                dd_weights[i][1] = -0.08333333333333333333;
                dd_weights[i][2] = +1.33333333333333333333;
                dd_weights[i][3] = -2.50000000000000000000;
                dd_weights[i][4] = +1.33333333333333333333;
                dd_weights[i][5] = -0.08333333333333333333;
            }
            else
            {
                dd_start[i] = -2;

                dd_weights[i][0] = -0.08333333333333333333;
                dd_weights[i][1] = +1.33333333333333333333;
                dd_weights[i][2] = -2.50000000000000000000;
                dd_weights[i][3] = +1.33333333333333333333;
                dd_weights[i][4] = -0.08333333333333333333;
                dd_weights[i][5] = 0;
            }
        }

    }

    out.du_stencil_start = d_start[0];
    out.dv_stencil_start = d_start[1];

    out.dudu_stencil_start = dd_start[0];
    out.dvdv_stencil_start = dd_start[1];
}
