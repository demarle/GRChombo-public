/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// code to calculate expansion jacobian for PETSc
// currently NOT being used
// instead, the code uses finite differences to calculate the derivative of the expansion
// using the code "AHFunction.inc"
// (TF): not sure if it works or not. Right now (10/2019), the 'form_jacobian' part of the
// algorithm takes only about 3% of the AH computation (using Chombo Timers).
// Whatever potential boost, seems uneeded
// Missing: converting all 'for (int a = 0; a < DIM; ++a)' for FOR1's FOR2's with CH_SPACEDIM

// =========================
//       XGEN :: BEGIN      
// =========================

// OUTPUTS
/*
double jac = 0;
double jacDu = 0;
double jacDv = 0;
double jacDuDu = 0;
double jacDvDv = 0;
double jacDuDv = 0;
*/

{

    // -------------------------
    //          Normal    
    // -------------------------
    
    Tensor<1, double> s = {0.};
    {
        for (int a = 0; a < DIM; ++a)
        {
            s[a] = data.df[a] - (deriv.duF * data.du[a]) - (deriv.dvF * data.dv[a]);
        }
    }
    
    Tensor<2, double> ds = {0.};
    {
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                ds[a][b] = data.ddf[a][b] - (deriv.duF * data.ddu[a][b]) - (deriv.dvF * data.ddv[a][b]) - (deriv.duduF * data.du[a] * data.du[b]) - (deriv.dvdvF * data.dv[a] * data.dv[b]) - (deriv.dudvF * (data.du[a] * data.dv[b] + data.du[b] * data.dv[a]));
            }
        }
    }
    
    Tensor<1, double> s_U = {0.};
    {
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                s_U[a] += data.g_UU[a][b] * s[b];
            }
        }
    }

    double normS = 0.0;
    {
        double _tmp0_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp0_ += data.g_UU[a][b] * s[a] * s[b];
            }
        }
        
        normS = sqrt(_tmp0_);
    }

    // -------------------------
    //    Jacobian of normal    
    // -------------------------

    Tensor<1, double> jacS = {0.};
    {
        for (int i = 0; i < DIM; ++i)
        {
            jacS[i] = data.ddfdf[i] - (data.ddudf[i] * deriv.duF) - (data.ddvdf[i] * deriv.dvF);
        }
    }
    
    Tensor<2, double> jacDs = {0.};
    {
        for (int i = 0; i < DIM; ++i)
        {
            for (int j = 0; j < DIM; ++j)
            {
                jacDs[i][j] = data.dddfdf[i][j] - (data.dddudf[i][j] * deriv.duF) - (data.dddvdf[i][j] * deriv.dvF) - (data.ddudf[j] * ((deriv.duduF * data.du[i]) + (deriv.dudvF * data.dv[i]))) - (data.ddvdf[j] * ((deriv.dudvF * data.du[i]) + (deriv.dvdvF * data.dv[i]))) - (((data.ddudf[i] * deriv.duduF) + (data.ddvdf[i] * deriv.dudvF)) * data.du[j]) - (((data.ddudf[i] * deriv.dudvF) + (data.ddvdf[i] * deriv.dvdvF)) * data.dv[j]);
            }
        }
    }

    Tensor<1, double> jacS_U = {0.};
    {
        for (int a = 0; a < DIM; ++a)
        {
            double _tmp0_ = 0.0;
            for (int b = 0; b < DIM; ++b)
            {
                _tmp0_ += data.g_UU[a][b] * jacS[b];
            }
            
            jacS_U[a] = _tmp0_;
        }
    }

    double tr_jacDs = {0.};
    {
        double _tmp0_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp0_ += data.g_UU[a][b] * jacDs[b][a];
            }
        }
        
        tr_jacDs = _tmp0_;
    }

    // -------------------------
    //       Extrinsic Aux    
    // -------------------------
    
    double tr_K = {0.};
    {
        double _tmp0_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp0_ += data.g_UU[a][b] * data.K[a][b];
            }
        }
        
        tr_K = _tmp0_;
    }
    
    Tensor<2, double> K_UU = {0.};
    {
        for (int b = 0; b < DIM; ++b)
        {
            for (int c = 0; c < DIM; ++c)
            {
                double _tmp0_ = 0.0;
                for (int a = 0; a < DIM; ++a)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        _tmp0_ += data.g_UU[a][b] * data.g_UU[c][d] * data.K[a][d];
                    }
                }
                
                K_UU[b][c] = _tmp0_;
            }
        }
    }
    
    Tensor<2, double> K_DU = {0.};
    {
        for (int b = 0; b < DIM; ++b)
        {
            for (int d = 0; d < DIM; ++d)
            {
                double _tmp0_ = 0.0;
                for (int a = 0; a < DIM; ++a)
                {
                    _tmp0_ += data.g_UU[a][d] * data.K[b][a];
                }
                
                K_DU[b][d] = _tmp0_;
            }
        }
    }
    
    // -------------------------
    //           Output      
    // -------------------------
    
    // Calculation of jac
    {
        double _tmp0_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp0_ += -2. * (normS * normS) * jacDs[a][b] * s_U[a] * s_U[b];
            }
        }
        
        double _tmp1_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            _tmp1_ += -2. * jacS_U[a] * tr_K * s[a];
        }
        
        double _tmp2_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    _tmp2_ += data.g_UU[b][c] * jacS_U[a] * data.dg[b][c][a];
                }
            }
        }
        
        double _tmp3_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    _tmp3_ += 2. * data.dxdf[a] * K_UU[b][c] * data.dg[b][c][a];
                }
            }
        }
        
        double _tmp4_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    _tmp4_ += -2. * data.dxdf[a] * data.g_UU[b][c] * data.dK[b][c][a];
                }
            }
        }
        
        double _tmp5_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    _tmp5_ += -2. * data.g_UU[b][c] * jacS_U[a] * data.dg[a][b][c];
                }
            }
        }
        
        double _tmp6_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int j = 0; j < DIM; ++j)
                    {
                        for (int k = 0; k < DIM; ++k)
                        {
                            _tmp6_ += -2. * data.dxdf[a] * data.g_UU[b][j] * data.g_UU[c][k] * data.dg[j][k][a] * ds[b][c];
                        }
                    }
                }
            }
        }
        
        double _tmp7_ = 0.0;
        for (int b = 0; b < DIM; ++b)
        {
            double _tmp8_ = 0.0;
            for (int c = 0; c < DIM; ++c)
            {
                double _tmp9_ = 0.0;
                for (int a = 0; a < DIM; ++a)
                {
                    _tmp9_ += -2. * jacS_U[a] * data.K[b][c] * s[a];
                }
                
                double _tmp10_ = 0.0;
                for (int a = 0; a < DIM; ++a)
                {
                    _tmp10_ += jacS_U[a] * data.dg[b][c][a];
                }
                
                double _tmp11_ = 0.0;
                for (int a = 0; a < DIM; ++a)
                {
                    _tmp11_ += data.dxdf[a] * tr_K * data.dg[b][c][a];
                }
                
                double _tmp12_ = 0.0;
                for (int a = 0; a < DIM; ++a)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        _tmp12_ += -4. * data.dxdf[a] * K_DU[b][d] * data.dg[c][d][a];
                    }
                }
                
                double _tmp13_ = 0.0;
                for (int a = 0; a < DIM; ++a)
                {
                    _tmp13_ += 2. * data.dxdf[a] * data.dK[b][c][a];
                }
                
                double _tmp14_ = 0.0;
                for (int a = 0; a < DIM; ++a)
                {
                    _tmp14_ += 2. * jacS_U[a] * data.dg[a][b][c];
                }
                
                double _tmp15_ = 0.0;
                for (int a = 0; a < DIM; ++a)
                {
                    _tmp15_ += 4. * jacS_U[a] * s[a] * ds[b][c];
                }
                
                double _tmp16_ = 0.0;
                for (int d = 0; d < DIM; ++d)
                {
                    double _tmp17_ = 0.0;
                    for (int a = 0; a < DIM; ++a)
                    {
                        for (int e = 0; e < DIM; ++e)
                        {
                            _tmp17_ += -2. * data.dxdf[a] * s_U[e] * data.dg[d][e][a] * ds[b][c];
                        }
                    }
                    
                    double _tmp18_ = 0.0;
                    for (int a = 0; a < DIM; ++a)
                    {
                        _tmp18_ += -2. * jacS_U[a] * s[a] * data.dg[b][c][d];
                    }
                    
                    double _tmp19_ = 0.0;
                    for (int a = 0; a < DIM; ++a)
                    {
                        for (int e = 0; e < DIM; ++e)
                        {
                            _tmp19_ += data.dxdf[a] * data.K[b][c] * s_U[e] * data.dg[d][e][a];
                        }
                    }
                    
                    double _tmp20_ = 0.0;
                    for (int a = 0; a < DIM; ++a)
                    {
                        for (int e = 0; e < DIM; ++e)
                        {
                            for (int i = 0; i < DIM; ++i)
                            {
                                _tmp20_ += -2. * data.dxdf[a] * data.g_UU[e][i] * data.dg[b][e][a] * data.dg[c][i][d];
                            }
                        }
                    }
                    
                    double _tmp21_ = 0.0;
                    for (int a = 0; a < DIM; ++a)
                    {
                        _tmp21_ += data.dxdf[a] * data.ddg[b][c][a][d];
                    }
                    
                    double _tmp22_ = 0.0;
                    for (int a = 0; a < DIM; ++a)
                    {
                        for (int e = 0; e < DIM; ++e)
                        {
                            for (int i = 0; i < DIM; ++i)
                            {
                                _tmp22_ += -1. * data.dxdf[a] * data.g_UU[e][i] * data.dg[d][i][a] * data.dg[b][c][e];
                            }
                        }
                    }
                    
                    double _tmp23_ = 0.0;
                    for (int a = 0; a < DIM; ++a)
                    {
                        for (int e = 0; e < DIM; ++e)
                        {
                            for (int i = 0; i < DIM; ++i)
                            {
                                _tmp23_ += data.dxdf[a] * s_U[e] * s_U[i] * data.dg[b][c][a] * data.dg[d][e][i];
                            }
                        }
                    }
                    
                    _tmp16_ += s_U[d] * (_tmp17_ + _tmp18_ + (normS * (_tmp19_ + (normS * (_tmp20_ + _tmp21_ + _tmp22_)))) + _tmp23_);
                }
                
                _tmp8_ += s_U[c] * ((normS * (_tmp9_ + (normS * (_tmp10_ + (normS * (_tmp11_ + _tmp12_ + _tmp13_)) + _tmp14_)))) + _tmp15_ + _tmp16_);
            }
            
            double _tmp24_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                _tmp24_ += -2. * jacS_U[a] * ds[b][a];
            }
            
            double _tmp25_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                _tmp25_ += -2. * jacS_U[a] * ds[a][b];
            }
            
            double _tmp26_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                _tmp26_ += 4. * jacS_U[a] * data.K[a][b];
            }
            
            double _tmp27_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        for (int e = 0; e < DIM; ++e)
                        {
                            for (int i = 0; i < DIM; ++i)
                            {
                                _tmp27_ += -1. * data.dxdf[a] * data.g_UU[c][d] * data.g_UU[e][i] * data.dg[c][e][a] * data.dg[d][i][b];
                            }
                        }
                    }
                }
            }
            
            double _tmp28_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        _tmp28_ += data.dxdf[a] * data.g_UU[c][d] * data.ddg[c][d][a][b];
                    }
                }
            }
            
            double _tmp29_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        for (int e = 0; e < DIM; ++e)
                        {
                            for (int i = 0; i < DIM; ++i)
                            {
                                _tmp29_ += -1. * data.dxdf[a] * data.g_UU[c][d] * data.g_UU[e][i] * data.dg[b][c][a] * data.dg[e][i][d];
                            }
                        }
                    }
                }
            }
            
            double _tmp30_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        _tmp30_ += -2. * data.dxdf[a] * data.g_UU[c][d] * data.ddg[b][c][a][d];
                    }
                }
            }
            
            double _tmp31_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        for (int e = 0; e < DIM; ++e)
                        {
                            for (int i = 0; i < DIM; ++i)
                            {
                                _tmp31_ += 2. * data.dxdf[a] * data.g_UU[c][d] * data.g_UU[e][i] * data.dg[d][i][a] * data.dg[b][c][e];
                            }
                        }
                    }
                }
            }
            
            double _tmp32_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        for (int e = 0; e < DIM; ++e)
                        {
                            for (int i = 0; i < DIM; ++i)
                            {
                                _tmp32_ += 2. * data.dxdf[a] * data.g_UU[c][d] * data.g_UU[e][i] * data.dg[b][c][a] * data.dg[d][e][i];
                            }
                        }
                    }
                }
            }
            
            double _tmp33_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        for (int j = 0; j < DIM; ++j)
                        {
                            for (int k = 0; k < DIM; ++k)
                            {
                                _tmp33_ += 2. * data.dxdf[a] * data.g_UU[c][j] * data.g_UU[d][k] * s[c] * data.dg[b][k][a] * ds[d][j];
                            }
                        }
                    }
                }
            }
            
            double _tmp34_ = 0.0;
            for (int a = 0; a < DIM; ++a)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        for (int m = 0; m < DIM; ++m)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                _tmp34_ += 2. * data.dxdf[a] * data.g_UU[c][d] * data.g_UU[m][n] * s[c] * data.dg[b][n][a] * ds[d][m];
                            }
                        }
                    }
                }
            }
            
            _tmp7_ += s_U[b] * (_tmp8_ + ((normS * normS) * (_tmp24_ + _tmp25_ + (normS * (_tmp26_ + (normS * (_tmp27_ + _tmp28_ + _tmp29_ + _tmp30_ + _tmp31_ + _tmp32_)))) + _tmp33_ + _tmp34_)));
        }
        
        jac = 0.5 * 1. / ((normS * normS) * (normS * normS)) * (_tmp0_ + ((normS * (normS * normS)) * (_tmp1_ + (normS * ((2. * tr_jacDs) + _tmp2_ + (normS * (_tmp3_ + _tmp4_)) + _tmp5_ + _tmp6_)))) + _tmp7_);
    }
    
    // Calculation of jacDu
    {
        double _tmp0_ = 0.0;
        for (int b = 0; b < DIM; ++b)
        {
            for (int c = 0; c < DIM; ++c)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    _tmp0_ += 2. * normS * data.K[c][d] * s_U[b] * s_U[c] * s_U[d] * data.du[b];
                }
            }
        }
        
        double _tmp1_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            double _tmp2_ = 0.0;
            for (int b = 0; b < DIM; ++b)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    _tmp2_ += -4. * (normS * (normS * normS)) * data.g_UU[b][d] * data.K[a][d] * data.du[b];
                }
            }
            
            double _tmp3_ = 0.0;
            for (int b = 0; b < DIM; ++b)
            {
                double _tmp4_ = 0.0;
                for (int c = 0; c < DIM; ++c)
                {
                    double _tmp5_ = 0.0;
                    for (int d = 0; d < DIM; ++d)
                    {
                        _tmp5_ += 2. * s_U[d] * data.du[a] * data.dg[b][c][d];
                    }
                    
                    _tmp4_ += s_U[c] * ((-4. * ds[a][b] * data.du[c]) + _tmp5_);
                }
                
                double _tmp6_ = 0.0;
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        _tmp6_ += -2. * data.g_UU[c][d] * data.dg[a][d][b] * data.du[c];
                    }
                }
                
                double _tmp7_ = 0.0;
                for (int c = 0; c < DIM; ++c)
                {
                    for (int e = 0; e < DIM; ++e)
                    {
                        _tmp7_ += -1. * data.g_UU[c][e] * data.du[c] * data.dg[a][b][e];
                    }
                }
                
                _tmp3_ += s_U[b] * (_tmp4_ + ((normS * normS) * ((2. * data.ddu[a][b]) + _tmp6_ + _tmp7_)));
            }
            
            _tmp1_ += s_U[a] * (_tmp2_ + _tmp3_);
        }
        
        double _tmp8_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp8_ += 2. * data.g_UU[a][b] * tr_K * s[a] * data.du[b];
            }
        }
        
        double _tmp9_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp9_ += -2. * data.g_UU[a][b] * data.ddu[a][b];
            }
        }
        
        double _tmp10_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        _tmp10_ += -1. * data.g_UU[a][c] * data.g_UU[b][d] * data.du[a] * data.dg[d][b][c];
                    }
                }
            }
        }
        
        double _tmp11_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int c = 0; c < DIM; ++c)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    for (int e = 0; e < DIM; ++e)
                    {
                        _tmp11_ += 2. * data.g_UU[a][c] * data.g_UU[d][e] * data.du[a] * data.dg[c][d][e];
                    }
                }
            }
        }
        
        double _tmp12_ = 0.0;
        for (int b = 0; b < DIM; ++b)
        {
            for (int c = 0; c < DIM; ++c)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    for (int e = 0; e < DIM; ++e)
                    {
                        _tmp12_ += 2. * data.g_UU[b][d] * data.g_UU[c][e] * s[b] * ds[c][d] * data.du[e];
                    }
                }
            }
        }
        
        double _tmp13_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int c = 0; c < DIM; ++c)
            {
                for (int i = 0; i < DIM; ++i)
                {
                    for (int j = 0; j < DIM; ++j)
                    {
                        _tmp13_ += 2. * data.g_UU[a][c] * data.g_UU[i][j] * s[a] * ds[c][i] * data.du[j];
                    }
                }
            }
        }
        
        jacDu = 0.5 * 1. / ((normS * normS) * (normS * normS)) * (_tmp0_ + _tmp1_ + ((normS * normS) * ((normS * (_tmp8_ + (normS * (_tmp9_ + _tmp10_ + _tmp11_)))) + _tmp12_ + _tmp13_)));
    }
    
    // Calculation of jacDv
    {
        double _tmp0_ = 0.0;
        for (int b = 0; b < DIM; ++b)
        {
            for (int c = 0; c < DIM; ++c)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    _tmp0_ += 2. * normS * data.K[c][d] * s_U[b] * s_U[c] * s_U[d] * data.dv[b];
                }
            }
        }
        
        double _tmp1_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            double _tmp2_ = 0.0;
            for (int b = 0; b < DIM; ++b)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    _tmp2_ += -4. * (normS * (normS * normS)) * data.g_UU[b][d] * data.K[a][d] * data.dv[b];
                }
            }
            
            double _tmp3_ = 0.0;
            for (int b = 0; b < DIM; ++b)
            {
                double _tmp4_ = 0.0;
                for (int c = 0; c < DIM; ++c)
                {
                    double _tmp5_ = 0.0;
                    for (int d = 0; d < DIM; ++d)
                    {
                        _tmp5_ += 2. * s_U[d] * data.dv[a] * data.dg[b][c][d];
                    }
                    
                    _tmp4_ += s_U[c] * ((-4. * ds[a][b] * data.dv[c]) + _tmp5_);
                }
                
                double _tmp6_ = 0.0;
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        _tmp6_ += -2. * data.g_UU[c][d] * data.dg[a][d][b] * data.dv[c];
                    }
                }
                
                double _tmp7_ = 0.0;
                for (int c = 0; c < DIM; ++c)
                {
                    for (int e = 0; e < DIM; ++e)
                    {
                        _tmp7_ += -1. * data.g_UU[c][e] * data.dv[c] * data.dg[a][b][e];
                    }
                }
                
                _tmp3_ += s_U[b] * (_tmp4_ + ((normS * normS) * ((2. * data.ddv[a][b]) + _tmp6_ + _tmp7_)));
            }
            
            _tmp1_ += s_U[a] * (_tmp2_ + _tmp3_);
        }
        
        double _tmp8_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp8_ += 2. * data.g_UU[a][b] * tr_K * s[a] * data.dv[b];
            }
        }
        
        double _tmp9_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp9_ += -2. * data.g_UU[a][b] * data.ddv[a][b];
            }
        }
        
        double _tmp10_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                for (int c = 0; c < DIM; ++c)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        _tmp10_ += -1. * data.g_UU[a][c] * data.g_UU[b][d] * data.dv[a] * data.dg[d][b][c];
                    }
                }
            }
        }
        
        double _tmp11_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int c = 0; c < DIM; ++c)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    for (int e = 0; e < DIM; ++e)
                    {
                        _tmp11_ += 2. * data.g_UU[a][c] * data.g_UU[d][e] * data.dv[a] * data.dg[c][d][e];
                    }
                }
            }
        }
        
        double _tmp12_ = 0.0;
        for (int b = 0; b < DIM; ++b)
        {
            for (int c = 0; c < DIM; ++c)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    for (int e = 0; e < DIM; ++e)
                    {
                        _tmp12_ += 2. * data.g_UU[b][d] * data.g_UU[c][e] * s[b] * ds[c][d] * data.dv[e];
                    }
                }
            }
        }
        
        double _tmp13_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int c = 0; c < DIM; ++c)
            {
                for (int i = 0; i < DIM; ++i)
                {
                    for (int j = 0; j < DIM; ++j)
                    {
                        _tmp13_ += 2. * data.g_UU[a][c] * data.g_UU[i][j] * s[a] * ds[c][i] * data.dv[j];
                    }
                }
            }
        }
        
        jacDv = 0.5 * 1. / ((normS * normS) * (normS * normS)) * (_tmp0_ + _tmp1_ + ((normS * normS) * ((normS * (_tmp8_ + (normS * (_tmp9_ + _tmp10_ + _tmp11_)))) + _tmp12_ + _tmp13_)));
    }
    
    // Calculation of jacDuDu
    {
        double _tmp0_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp0_ += -1. * (normS * normS) * data.g_UU[a][b] * data.du[a] * data.du[b];
            }
        }
        
        double _tmp1_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp1_ += s_U[a] * s_U[b] * data.du[a] * data.du[b];
            }
        }
        
        jacDuDu = 1. / (normS * normS) * (_tmp0_ + _tmp1_);
    }
    
    // Calculation of jacDvDv
    {
        double _tmp0_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp0_ += -1. * (normS * normS) * data.g_UU[a][b] * data.dv[a] * data.dv[b];
            }
        }
        
        double _tmp1_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp1_ += s_U[a] * s_U[b] * data.dv[a] * data.dv[b];
            }
        }
        
        jacDvDv = 1. / (normS * normS) * (_tmp0_ + _tmp1_);
    }
    
    // Calculation of jacDuDv
    {
        double _tmp0_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp0_ += -2. * (normS * normS) * data.g_UU[a][b] * data.dv[a] * data.du[b];
            }
        }
        
        double _tmp1_ = 0.0;
        for (int a = 0; a < DIM; ++a)
        {
            for (int b = 0; b < DIM; ++b)
            {
                _tmp1_ += 2. * s_U[a] * s_U[b] * data.du[a] * data.dv[b];
            }
        }
        
        jacDuDv = 1. / (normS * normS) * (_tmp0_ + _tmp1_);
    }
}

// =========================
//        XGEN :: END       
// =========================
