/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// code included in ApparentHorizon.impl.hpp::FormFunction and FormJacobian to calculate expansion

// =========================
//       XGEN :: BEGIN      
// =========================

// OUTPUTS
/*
double expansion = 0;
*/

{
    
    Tensor<1, double> s = {0.};
    FOR1(a) { s[a] = data.df[a] - (deriv.duF * data.du[a]) - (deriv.dvF * data.dv[a]); }
    
    Tensor<2, double> ds = {0.};
    FOR2(a,b)
    {
        ds[a][b] = data.ddf[a][b] - (deriv.duF * data.ddu[a][b]) - (deriv.dvF * data.ddv[a][b]) - (deriv.duduF * data.du[a] * data.du[b]) - (deriv.dvdvF * data.dv[a] * data.dv[b]) - (deriv.dudvF * (data.du[a] * data.dv[b] + data.du[b] * data.dv[a]));
    }
    
    Tensor<1, double> s_U = {0.};
    FOR2(a,b) { s_U[a] += data.g_UU[a][b] * s[b]; }

    double norm_s = 0.0;
    FOR2(a,b) { norm_s += data.g_UU[a][b] * s[a] * s[b]; }
    norm_s = sqrt(norm_s);

    Tensor<3, double> chris = {0.};
    FOR4(a,b,c,d)
    {
        chris[a][b][c] += 0.5 * data.g_UU[a][d] * (data.dg[b][d][c] + data.dg[c][d][b] - data.dg[b][c][d]);
    }

    Tensor<1, double> S = {0.};
    Tensor<1, double> S_U = {0.};
    Tensor<2, double> DS = {0.};
    FOR1(a)
    {
        S[a] = s[a] / norm_s;
        S_U[a] = s_U[a] / norm_s;

        FOR1(b)
        {
            DS[a][b] = ds[a][b] / norm_s;

            FOR1(c)
            {
                DS[a][b] -= (chris[c][a][b] * s[c] + ds[a][c] * s[b] * s_U[c] / (norm_s * norm_s)) / norm_s;

                FOR1(d)
                {
                    DS[a][b] += (chris[d][a][c] * s[b] * s_U[c] * s[d]) / (norm_s * norm_s * norm_s);
                }
            }
        }
    }
    
    // -------------------------
    //          Output      
    // -------------------------
    
    // Calculation of expansion
    expansion = 0;
    FOR2(a,b) { expansion += (data.g_UU[a][b] - S_U[a] * S_U[b]) * (DS[a][b] - data.K[a][b]); }
}

// =========================
//        XGEN :: END       
// =========================
