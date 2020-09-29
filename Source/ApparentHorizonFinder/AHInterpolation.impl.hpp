/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHINTERPOLATION_IMPL_HPP_
#define _AHINTERPOLATION_IMPL_HPP_

#if !defined(_AHINTERPOLATION_HPP_)
#error "This file should only be included through AHInterpolation.hpp"
#endif

template <typename CoordSystem>
AHInterpolation<CoordSystem>::AHInterpolation(
    const CoordSystem &a_coord_system,
    AMRInterpolator<Lagrange<4>> *a_interpolator)
    : m_coord_system(a_coord_system), m_interpolator(a_interpolator)
{
    ;

    // below:
    // determine maximum and minimum physical coordinate of the grid
    // (so that 'fit_in_grid' knows when PETSc has diverged out of the grid and
    // doesn't let him do so,
    //  as it would cause an error in the AMRInterpolator)

    const Box &domainBox = const_cast<AMR &>(m_interpolator->getAMR())
                               .getAMRLevels()[0]
                               ->problemDomain()
                               .domainBox();
    const IntVect &small_end = domainBox.smallEnd();
    const IntVect &big_end = domainBox.bigEnd();

    const std::array<double, CH_SPACEDIM> &coarsest_dx =
        m_interpolator->get_coarsest_dx();
    const std::array<double, CH_SPACEDIM> &coarsest_origin =
        m_interpolator->get_coarsest_origin();

    // set coordinate minimum and maximum as ONE cells before the boundary
    for (unsigned i = 0; i < CH_SPACEDIM; ++i)
    {
        if (m_interpolator->get_boundary_reflective(Side::Lo, i))
        {
            m_coord_min[i] =
                -(big_end[i] - 1) * coarsest_dx[i] - coarsest_origin[i];
        }
        else
        {
            m_coord_min[i] =
                (small_end[i] + 1) * coarsest_dx[i] + coarsest_origin[i];
        }

        if (m_interpolator->get_boundary_reflective(Side::Hi, i))
        {
            m_coord_max[i] =
                (2 * big_end[i] - 1) * coarsest_dx[i] + coarsest_origin[i];
        }
        else
        {
            m_coord_max[i] =
                (big_end[i] - 1) * coarsest_dx[i] + coarsest_origin[i];
        }
        // pout() << "m_coord = " << m_coord_min[i] << "\t" << m_coord_max[i]
        //       << " with coarse origin " << coarsest_origin[i] << " and dx "
        //       << coarsest_dx[i] << std::endl;
    }
}

template <typename CoordSystem>
void AHInterpolation<CoordSystem>::set_SNES(SNES &a_snes)
{
    m_snes = &a_snes;
}

template <typename CoordSystem>
const AMRInterpolator<Lagrange<4>> *
AHInterpolation<CoordSystem>::get_interpolator() const
{
    return m_interpolator;
}

template <typename CoordSystem>
array<bool, CH_SPACEDIM - 1> AHInterpolation<CoordSystem>::is_periodic() const
{
    return m_coord_system.is_periodic();
}

template <typename CoordSystem>
array<double, CH_SPACEDIM - 1>
AHInterpolation<CoordSystem>::get_domain_low() const
{
    return m_coord_system.get_domain_low();
}
template <typename CoordSystem>
array<double, CH_SPACEDIM - 1>
AHInterpolation<CoordSystem>::get_domain_high() const
{
    return m_coord_system.get_domain_high();
}

template <typename CoordSystem>
std::vector<std::string> AHInterpolation<CoordSystem>::get_labels() const
{
    return m_coord_system.get_labels();
}

template <typename CoordSystem>
void AHInterpolation<CoordSystem>::set_center(
    const std::array<double, CH_SPACEDIM> &center)
{
    m_coord_system.set_center(center);
}

template <typename CoordSystem>
bool AHInterpolation<CoordSystem>::fit_in_grid(double &x, double &y, double &z)
{
    CH_TIME("AHInterpolation::fit_in_grid");

    bool out_of_grid = false;

    // if out of bounds, put back in grid
    if (x < m_coord_min[0])
    {
        out_of_grid = true;
        x = m_coord_min[0];
    }
    if (x > m_coord_max[0])
    {
        out_of_grid = true;
        x = m_coord_max[0];
    }

    if (y < m_coord_min[1])
    {
        out_of_grid = true;
        y = m_coord_min[1];
    }
    if (y > m_coord_max[1])
    {
        out_of_grid = true;
        y = m_coord_max[1];
    }

    if (z < m_coord_min[2])
    {
        out_of_grid = true;
        z = m_coord_min[2];
    }
    if (z > m_coord_max[2])
    {
        out_of_grid = true;
        z = m_coord_max[2];
    }

    return out_of_grid;
}

template <typename CoordSystem>
array<double, 3> AHInterpolation<CoordSystem>::transform(double u, double v,
                                                         double f)
{
    CH_TIME("AHInterpolation::resize_vectors");

    return m_coord_system.transform(u, v, f);
}

template <typename CoordSystem>
void AHInterpolation<CoordSystem>::resize_vectors(int n)
{
    CH_TIME("AHInterpolation::resize_vectors");

    m_u.resize(n);
    m_v.resize(n);
    m_f.resize(n);

    m_x.resize(n);
    m_y.resize(n);
    m_z.resize(n);

    for (int i = 0; i < NUM_VARS; ++i)
    {
        m_fun[i].resize(n);
    }

    for (int i = c_chi; i <= c_h33; ++i)
    {
        for (int j = 0; j < CH_SPACEDIM; ++j)
        {
            m_dfun[j][i].resize(n);
        }
    }
}

//! triplet of functions to be used together in blocks of code that require
//! PETSc AND AMRInterpolator to interpolate 'keep_interpolating_if_inactive'
//! returns 'true' immediately for PETSc ranks for non-PETSc ranks, it loops in
//! a sequence of 'interpolate()' waiting for PETSc ranks to call
//! 'interpolate()'' as well (as this needs to be called by all Chombo
//! processes) can be aborted by calling 'break_interpolation_loop()' for PETSc
//! ranks
/** Example of usage:
 * if(m_geom.keep_interpolating_if_inactive())
 * {
 *     (... PETSc code that calls 'set_data()'...)
 *     m_geom.break_interpolation_loop();
 * }
 */
template <typename CoordSystem>
bool AHInterpolation<CoordSystem>::keep_interpolating_if_inactive()
{
    CH_TIME("AMRInterpolator::keep_interpolating_if_inactive");

    if (!AHFinder::is_rank_active())
    {
        int keep_interpolating = 1;
        while (keep_interpolating)
            keep_interpolating = interpolate();
        return false;
    }
    else
        return true;
}

template <typename CoordSystem>
void AHInterpolation<CoordSystem>::break_interpolation_loop() const
{
    CH_TIME("AMRInterpolator::break_interpolation_loop");

    // break "keep_interpolating" loop
    int ZERO = 0;
    int keep_interpolating = 0;
#ifdef CH_MPI
    // this is also called in 'interpolate()', breaking the loop
    MPI_Allreduce(&ZERO, &keep_interpolating, 1, MPI_INT, MPI_LAND,
                  Chombo_MPI::comm);
#endif
}

template <typename CoordSystem> int AHInterpolation<CoordSystem>::interpolate()
{
    CH_TIME("AHInterpolation::interpolate");

    // Code below used to allow PETSc to run in a sub-communicator.
    // For that, PETSc processes run 'set_data' and 'interpolate', while
    // non-PETSc processes loop through 'interpolate'. If all cores run
    // interpolate, 'MPI_Allreduce' will return 'keep_interpolating' as true. To
    // exit the loop for non-PETSc cores, simply do on them a call to:
    // MPI_Allreduce(&ZERO, &keep_interpolating, 1, MPI_INT, MPI_LAND,
    // Chombo_MPI::comm); (this is done in 'break_interpolation_loop()')

    int ONE = 1;
    int keep_interpolating = 1;
#ifdef CH_MPI
    MPI_Allreduce(&ONE, &keep_interpolating, 1, MPI_INT, MPI_LAND,
                  Chombo_MPI::comm);
#endif

    if (!keep_interpolating)
        return 0;

    InterpolationQuery query(m_x.size());
    query.setCoords(0, &m_x[0]).setCoords(1, &m_y[0]).setCoords(2, &m_z[0]);

    for (int i = 0; i < NUM_VARS; ++i)
    {
        query.addComp(i, &m_fun[i][0]);
    }

    // Derivatives required only for the metric component
    for (int i = c_chi; i <= c_h33; ++i)
    {
        query.addComp(i, &m_dfun[0][i][0], Derivative::dx)
            .addComp(i, &m_dfun[1][i][0], Derivative::dy)
            .addComp(i, &m_dfun[2][i][0], Derivative::dz);
    }

    m_interpolator->interp(query);

    return 1;
}

template <typename CoordSystem>
void AHInterpolation<CoordSystem>::refresh_interpolator()
{
    CH_TIME("AHInterpolation::refresh_interpolator");
    m_interpolator->refresh();
}

// 'set_data' calls 'interpolate'. All Chombo_MPI:comm need to run 'interpolate'
// for it to work
template <typename CoordSystem>
void AHInterpolation<CoordSystem>::set_data(const vector<double> &u,
                                            const vector<double> &v,
                                            const vector<double> &f)
{
    CH_TIME("AHInterpolation::set_data");

    // TODO: The following code only work when CH_SPACEDIM = 3;
    CH_assert(CH_SPACEDIM == 3);

    CH_assert(u.size() == v.size());
    CH_assert(u.size() == f.size());

    const int n = f.size();
    resize_vectors(n);

    bool out_of_grid = false;

    // Transform to Cartesian
    for (int i = 0; i < n; ++i)
    {
        m_u[i] = u[i];
        m_v[i] = v[i];
        m_f[i] = f[i];

        array<double, 3> xyz = m_coord_system.transform(u[i], v[i], f[i]);
        m_x[i] = xyz[0];
        m_y[i] = xyz[1];
        m_z[i] = xyz[2];

        // don't let PETSc diverge to outside of the grid (this can happen if
        // there is no BH)
        out_of_grid |= fit_in_grid(m_x[i], m_y[i], m_z[i]);
    }

    interpolate();

    // abort if out of grid - reduces the time to diverge dramatically
    if (out_of_grid)
        SNESSetFunctionDomainError(*m_snes);
}

template <typename CoordSystem>
AHGeometryData AHInterpolation<CoordSystem>::get_data(int idx) const
{
    CH_TIME("AHInterpolation::get_data");

    // TODO: The following code only work when CH_SPACEDIM = 3;
    CH_assert(CH_SPACEDIM == 3);

    AHGeometryData out;
    out.lapse = m_fun[c_lapse][idx];

    // * ---------------------------
    // *   COORDINATE-RELATED DATA
    // * ---------------------------

    m_coord_system.transform_derivs(m_u[idx], m_v[idx], m_f[idx], out);

    // * ---------------------------
    // *       GR-RELATED DATA
    // * ---------------------------

    const int h[3][3] = {
        {c_h11, c_h12, c_h13}, {c_h12, c_h22, c_h23}, {c_h13, c_h23, c_h33}};

    const int A[3][3] = {
        {c_A11, c_A12, c_A13}, {c_A12, c_A22, c_A23}, {c_A13, c_A23, c_A33}};

    const double chi = m_fun[c_chi][idx];
    const double trK = m_fun[c_K][idx];

    // INVERSE METRIC
    Tensor<2, double, CH_SPACEDIM> h_DD;
    FOR2(i, j) { h_DD[i][j] = m_fun[h[i][j]][idx]; }

    Tensor<2, double, CH_SPACEDIM> h_UU =
        TensorAlgebra::compute_inverse_sym(h_DD);
    FOR2(i, j) { out.g_UU[i][j] = chi * h_UU[i][j]; }

    // Reconstructing ADM variables

    for (int i = 0; i < CH_SPACEDIM; ++i)
    {
        for (int j = i; j < CH_SPACEDIM; ++j)
        {
            {
                const double g = m_fun[h[i][j]][idx] / (chi);
                out.g[i][j] = g;
                out.g[j][i] = g;
            }

            {
                const double K =
                    m_fun[A[i][j]][idx] / chi + trK * out.g[i][j] / CH_SPACEDIM;
                out.K[i][j] = K;
                out.K[j][i] = K;
            }

            for (int k = 0; k < CH_SPACEDIM; ++k)
            {
                {
                    const double dg =
                        (m_dfun[k][h[i][j]][idx] -
                         (m_fun[h[i][j]][idx] * m_dfun[k][c_chi][idx]) / chi) /
                        chi;
                    out.dg[i][j][k] = dg;
                    out.dg[j][i][k] = dg;
                }
            }
        }
    }

    out.trK = trK;
    return out;
}

#endif /* _AHINTERPOLATION_IMPL_HPP_ */
