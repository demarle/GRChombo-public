/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// code included in ApparentHorizon.impl.hpp::CalculateArea to calculate area element

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

    
    Tensor<1, double> s_U = {0.};
    FOR2(a,b) { s_U[a] += data.g_UU[a][b] * s[b]; }

    double norm_s = 0.0;
    FOR2(a,b) { norm_s += data.g_UU[a][b] * s[a] * s[b]; }
    norm_s = sqrt(norm_s);


    Tensor<1, double> S = {0.};
    Tensor<1, double> S_U = {0.};
    FOR1(a)
    {
        S[a] = s[a] / norm_s;
        S_U[a] = s_U[a] / norm_s;
    }
    
    // -------------------------
    //          Output      
    // -------------------------
    
    // Calculation of expansion
    {
        Tensor<2, double> mat = {0.};
	    double lambda1 = 0;
	    double lambda2 = 0;
	    double root = 0;
	
        FOR2(a,b) { mat[a][b] = (data.g_UU[a][b] - S_U[a] * S_U[b]); }


	// Getting Eigenvalues for matrix with det(mat) = 0
	root = mat[0][0] * mat[0][0] + 4 * mat[0][1] * mat[1][0]
           - 2 * mat[0][0] * mat[1][1] + mat[1][1] * mat[1][1]
           + 4 * mat[0][2] * mat[2][0] + 4 * mat[1][2] * mat[2][1]
           - 2 * mat[0][0] * mat[2][2] - 2 * mat[1][1] * mat[2][2]
           + mat[2][2] * mat[2][2];
	// Cutoff for small values, to avoid NaNs from sqrt
	if (abs(root) < 1e-16 || root < 0.0 ){
		root = 0.0;
	}
    
    //pout() << tag << "root = " << root << "\n";
    lambda1 = (mat[0][0] + mat[1][1] + mat[2][2])/2. - sqrt(root)/2.;
	lambda2 = (mat[0][0] + mat[1][1] + mat[2][2])/2. + sqrt(root)/2.;
    /*
    pout() << tag << "lambda1 = " << lambda1 << ", lambda2 = " << lambda2 
           << "\n";
    */
	area_element = 1.0/(lambda1*lambda2);

   }
}

// =========================
//        XGEN :: END       
// =========================
