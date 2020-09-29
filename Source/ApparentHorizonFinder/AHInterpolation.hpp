/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHINTERPOLATION_HPP_
#define _AHINTERPOLATION_HPP_

#include "AHFinder.hpp"
#include "AHGeometryData.hpp"
#include "AMRInterpolator.hpp"
#include "Lagrange.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"

//! Class used for interpolation of the variables needed to calculate expansion
//! with the data from a given 'CoordSystem'
template <typename CoordSystem> class AHInterpolation
{
  private:
    CoordSystem m_coord_system;
    AMRInterpolator<Lagrange<4>> *m_interpolator;

    // variables of AH in 'CoordSystem'
    // (in Spherical Coordinates Class 'AHSphericalCoords',
    // they correspond to 'theta', 'phi' and 'sqrt(radius)')
    vector<double> m_u;
    vector<double> m_v;
    vector<double> m_f;

    // variables of AH in cartesian (for AMR interpolation)
    vector<double> m_x;
    vector<double> m_y;
    vector<double> m_z;

    //! vector with GRChombo variables (chi, h, K)
    vector<double> m_fun[NUM_VARS];
    //! vector with derivatives of chi and 3-metric
    vector<double> m_dfun[CH_SPACEDIM][c_h33 - c_chi + 1];

    void resize_vectors(int n);

    std::array<double, CH_SPACEDIM> m_coord_min,
        m_coord_max; //!< maximum and minimum of level 0 box, used in
                     //!< 'fit_in_grid'

    SNES *m_snes;

    //! when PETSc tried to diverge out of the grid, this doesn't let him do so
    //! it forces him to stay on the grid, cause non-convergence
    bool fit_in_grid(double &x, double &y, double &z);

  public:
    AHInterpolation(const CoordSystem &a_coordSystem,
                    AMRInterpolator<Lagrange<4>> *a_interpolator);

    void set_SNES(SNES &a_snes);

    const AMRInterpolator<Lagrange<4>> *get_interpolator() const;

    // several of the methods below just call the correspondent method of the
    // CoordSystem, as these are needed in the 'ApparentHorizon' class
    array<bool, CH_SPACEDIM - 1> is_periodic() const;
    array<double, CH_SPACEDIM - 1>
    get_domain_low() const; //!< lower bound of (u,v) coordinates
    array<double, CH_SPACEDIM - 1>
    get_domain_high() const; //!< upper bound of (u,v) coordinates
    std::vector<std::string>
    get_labels() const; //!< transform from CoordSystem to Cartesian
    void set_center(
        const std::array<double, CH_SPACEDIM> &); //!< set center of CoordSystem
    void refresh_interpolator(); //!< refresh AMRInterpolator 'm_interpolator'

    array<double, 3>
    transform(double u, double v,
              double f); //!< transform from spherical to cartesian
    array<double, 3>
    transform_inverse(double u, double v,
                      double f); //!< transform derivatives from spherical
                                 //!< to cartesian derivatives

    //! 'set_data' calls 'interpolate'. All Chombo_MPI:comm need to run
    //! 'interpolate' for it to work (if some are not part of PETSc, they still
    //! need to run 'interpolate' on the side) see further description below
    void set_data(const vector<double> &u, const vector<double> &v,
                  const vector<double> &f);
    AHGeometryData get_data(int idx) const;

    //! triplet of functions to be used together in blocks of code that require
    //! PETSc AND AMRInterpolator to interpolate
    //! 'keep_interpolating_if_inactive' returns 'true' immediately for PETSc
    //! ranks for non-PETSc ranks, it loops in a sequence of 'interpolate()'
    //! waiting for PETSc ranks to call 'interpolate()'' as well (as this needs
    //! to be called by all Chombo processes) can be aborted by calling
    //! 'break_interpolation_loop()' for PETSc ranks
    /** Example of usage:
     * if(m_geom.keep_interpolating_if_inactive())
     * {
     *     (... PETSc code that calls 'set_data()'...)
     *     m_geom.break_interpolation_loop();
     * }
     */
    bool keep_interpolating_if_inactive();
    void
    break_interpolation_loop() const; //!< see 'keep_interpolating_if_inactive'
                                      //!< (this one could be static)
    int interpolate();                //!< see 'keep_interpolating_if_inactive'
};

#include "AHInterpolation.impl.hpp"

#endif /* _AHINTERPOLATION_HPP_ */
