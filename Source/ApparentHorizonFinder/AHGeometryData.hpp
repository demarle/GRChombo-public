/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHGEOMETRYDATA_HPP_
#define _AHGEOMETRYDATA_HPP_

#include "Tensor.hpp"

// The d prefix refers to partial derivatives wrt Cartesian coordinates

struct AHGeometryData
{

    // transform jacobian
    Tensor<1, double> du;
    Tensor<1, double> dv;
    Tensor<1, double> df;

    // transform hessian
    Tensor<2, double> ddu;
    Tensor<2, double> ddv;
    Tensor<2, double> ddf;

    // inverse jacobian
    Tensor<1, double> dxdu;
    Tensor<1, double> dxdv;
    Tensor<1, double> dxdf;

    // inverse hessian

    // induced metric
    Tensor<2, double> g;
    Tensor<2, double> g_UU;
    Tensor<3, double> dg;

    // extrinsic curvature
    Tensor<2, double> K;

    double trK;
    double lapse;

    AHGeometryData()
    {
        // force all (double) elements of AHGeometryData to be 0
        memset(this, 0, sizeof(AHGeometryData));
    }
};

#endif /* _AHGEOMETRYDATA_HPP_ */
