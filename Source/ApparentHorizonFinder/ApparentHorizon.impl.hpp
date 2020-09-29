/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _APPARENTHORIZON_IMPL_HPP_
#define _APPARENTHORIZONI_IMPL_HPP_

#if !defined(_APPARENTHORIZON_HPP_)
#error "This file should only be included through ApparentHorizon.hpp"
#endif

#include "AHGeometryData.hpp"
#include "AMR.H"
#include "SmallDataIO.hpp"
#include "TensorAlgebra.hpp"
#include "petscsys.h"
#include "petscviewerhdf5.h"
#include <algorithm> // for std::max_element
#include <iomanip>   // for 'setw' in stringstream

template <typename AHGeom>
ApparentHorizon<AHGeom>::ApparentHorizon(
    AHGeom &a_geom, const std::array<double, CH_SPACEDIM> &a_center,
    double a_initial_guess, const AHFinder::params &a_params,
    const std::string &a_stats, const std::string &a_coords)
    : m_geom(a_geom), m_geom_plus(a_geom), m_geom_minus(a_geom),

      m_params(a_params),

      m_stats(a_stats), m_coords(a_coords),

      m_converged(false), m_has_been_found(false),

      m_num_failed_convergences(0),

      m_center(a_center), // set below in 'set_center(center)'

      m_initial_guess(a_initial_guess),

      m_merger_pair({-1, -1}),

      m_periodic_u(a_geom.is_periodic()[0]),

      m_num_global_u(a_params.num_points_u),
      m_num_global_v(a_params.num_points_v),

#if CH_SPACEDIM == 3
      m_periodic_v(m_geom.is_periodic()[1])
#elif CH_SPACEDIM == 2
      m_periodic_v(false)
#endif
{
    initialise_PETSc();
    set_center(a_center);

    int current_step = AMR::s_step * pow(2.0, m_params.level);

    if (current_step == 0)
        // solve for first time step if t==0
        solve(1., 0., 0.,
              true); // 'dt' doesn't matter for 1st time step, so set to 1
    else
    {
        double time = m_geom.get_interpolator()->getAMR().getCurrentTime();
        restart(current_step, time);
        // double dt = time / (double)current_step;
        // write_coords_to_file(dt, time, time);
    }
}

template <typename AHGeom> ApparentHorizon<AHGeom>::~ApparentHorizon()
{
    finalise_PETSc();
}

template <typename AHGeom>
bool ApparentHorizon<AHGeom>::do_solve(double a_dt, double a_time) const
{
    CH_assert(a_dt != 0); // Check if time was set!
    return !(((int)(a_time / a_dt)) % m_params.solve_interval);
}
template <typename AHGeom>
bool ApparentHorizon<AHGeom>::do_print(double a_dt, double a_time) const
{
    CH_assert(a_dt != 0); // Check if time was set!
    return !(((int)(a_time / a_dt)) %
             (m_params.solve_interval * m_params.print_interval));
}

template <typename AHGeom>
std::array<double, CH_SPACEDIM> ApparentHorizon<AHGeom>::get_center() const
{
    return m_center;
}

template <typename AHGeom> double ApparentHorizon<AHGeom>::get_max_F() const
{
    return m_max_F;
}

template <typename AHGeom> double ApparentHorizon<AHGeom>::get_min_F() const
{
    return m_min_F;
}

template <typename AHGeom>
void ApparentHorizon<AHGeom>::set_center(
    const std::array<double, CH_SPACEDIM> &a_center)
{
    m_center_old = m_center;
    m_center = a_center;
    m_geom.set_center(m_center);
    m_geom_plus.set_center(m_center);
    m_geom_minus.set_center(m_center);
}

template <typename AHGeom> bool ApparentHorizon<AHGeom>::get_converged() const
{
    return m_converged;
}

template <typename AHGeom> bool ApparentHorizon<AHGeom>::has_been_found() const
{
    return m_has_been_found;
}

template <typename AHGeom> bool ApparentHorizon<AHGeom>::stop_solving() const
{
    return (m_has_been_found &&
            m_num_failed_convergences > m_params.max_failed_convergences &&
            m_params.max_failed_convergences != 0);
}

template <typename AHGeom> double ApparentHorizon<AHGeom>::get_initial_guess()
{
    return m_initial_guess;
}

template <typename AHGeom>
void ApparentHorizon<AHGeom>::set_merger_pair(int first, int second)
{
    CH_assert(first >= 0 && second >= 0);
    m_merger_pair.first = first;
    m_merger_pair.second = second;
}
template <typename AHGeom>
std::pair<int, int> ApparentHorizon<AHGeom>::get_merger_pair()
{
    return m_merger_pair;
}
template <typename AHGeom> bool ApparentHorizon<AHGeom>::is_merger()
{
    return m_merger_pair.first >= 0; // if ==-1, not merger
}

template <typename AHGeom> void ApparentHorizon<AHGeom>::set_initial_guess()
{
    CH_TIME("ApparentHorizon::set_initial_guess");

    if (!AHFinder::is_rank_active())
        return;

    pout() << "Setting Initial Guess to f=" << m_initial_guess << std::endl;

    // read PETSc array to 'f'
    dmda_arr_t f;
    DMDAVecGetArray(m_dmda, m_snes_soln, &f);

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
#endif
    {
        for (int u = m_umin; u < m_umax; ++u)
        {
#if CH_SPACEDIM == 3
            double &f_point = f[v][u];
#else
            double &f_point = f[u];
#endif

            f_point = sqrt(
                m_initial_guess); // f = sqrt(radius) for spherical coordinates
        }
    }

    // write PETSc array back
    DMDAVecRestoreArray(m_dmda, m_snes_soln, &f);
}

template <typename AHGeom>
void ApparentHorizon<AHGeom>::solve(double a_dt, double a_time,
                                    double a_restart_time, bool a_first_step)
{
    CH_TIME("ApparentHorizon::solve");

    if (!do_solve(a_dt, a_time))
        return;

    m_geom.refresh_interpolator(); //(ALL CHOMBO ranks do it!!)

    // PETSc processes go inside 'if', others "wait" until 'if' gets to
    // 'm_geom.break_interpolation_loop()'
    if (m_geom.keep_interpolating_if_inactive())
    {
        CH_TIME("ApparentHorizon::solve::solving");

        if (!get_converged())
            set_initial_guess(); // reset initial guess if diverged (or in first
                                 // attempt)

        // actual solve happens here!
        SNESSolve(m_snes, NULL, m_snes_soln);

        PetscInt its;
        SNESGetIterationNumber(m_snes, &its);
        pout() << "SNES Iteration number " << its << endl;
        SNESGetLinearSolveIterations(m_snes, &its);
        pout() << "KSP Iteration number " << its << endl;

        m_geom.break_interpolation_loop();
    }

    double area, spin, mass, coord_area;
    {
        CH_TIME("ApparentHorizon::solve::post-solving");

        // ask PETSc if it converged
        check_convergence();

        area = calculate_area();
        spin = calculate_spin(area);
        mass = sqrt(area / (8.0 * M_PI * (1 + sqrt(1 - spin * spin))));
        coord_area = calculate_coord_area();
        if (m_params.updateCenter)
            calculate_center();
        if (m_converged)
        {
            calculate_minmax_F();
        }
        else
        {
            m_max_F = 0.0;
            m_min_F = 0.0;
        }
    }

    // print stats (area, spin, center) and coordinates
    // stop printing if it stopped converging
    if ((m_converged || !m_has_been_found) && do_print(a_dt, a_time))
    {
        CH_TIME("ApparentHorizon::solve::printing");
        pout() << "Printing statistics and coordinates." << std::endl;

        int pre = 4;
        std::vector<double> values(pre + CH_SPACEDIM);

        values[0] = area;
        values[1] = spin;
        values[2] = mass;
        values[3] = coord_area;
        for (int i = 0; i < CH_SPACEDIM; ++i)
            values[pre + i] = m_center[i];

        // write stats

        // real_dt is relevant for the 'remove_duplicate_time_data' part
        double real_dt =
            a_dt * m_params.solve_interval * m_params.print_interval;
        SmallDataIO file(m_stats, real_dt, a_time, a_restart_time,
                         SmallDataIO::APPEND, a_first_step);

        file.remove_duplicate_time_data();

        // print headers to stats file in the beginning of evolution
        if (a_time == 0.0)
        {
            std::vector<std::string> headers(pre + CH_SPACEDIM);

            headers[0] = "area";
            headers[1] = "spin";
            headers[2] = "mass";
            headers[3] = "coord area";
            headers[pre] = "c_x";
            headers[pre + 1] = "c_y";
            headers[pre + 2] = "c_z";

            file.write_header_line(headers);
        }

        file.write_time_data_line(values);

        // write coordinates
        write_coords_to_file(a_dt, a_time, a_restart_time);
    }
}

template <typename AHGeom> void ApparentHorizon<AHGeom>::check_convergence()
{
    CH_TIME("ApparentHorizon::check_convergence");

    int result;
    if (AHFinder::is_rank_active())
    {
        SNESConvergedReason reason;
        SNESGetConvergedReason(m_snes, &reason);

        // result will be 0 if any of the PETSc ranks says 'reason <=0' (PETSc
        // convergence error)
        int test = (reason > 0);
#ifdef CH_MPI
        MPI_Allreduce(&test, &result, 1, MPI_INT, MPI_LAND, Chombo_MPI::comm);
#else
        result = test;
#endif
    }
    else
    {
        int ONE = 1;
#ifdef CH_MPI
        MPI_Allreduce(&ONE, &result, 1, MPI_INT, MPI_LAND, Chombo_MPI::comm);
#else
        result = ONE;
#endif
    }

    m_converged = (bool)result;
    if (m_converged)
    {
        m_num_failed_convergences = 0;
    }
    else
    {
        ++m_num_failed_convergences;
    }

    pout() << (m_converged ? "Solver converged. Horizon FOUND."
                           : "Solver diverged. Horizon NOT found.")
           << std::endl;

    if (!m_has_been_found && m_converged)
        m_has_been_found = true; // finally found :D
}

template <typename AHGeom>
void ApparentHorizon<AHGeom>::restart(int a_int_step, double a_current_time)
{
    CH_TIME("ApparentHorizon::restart");

    // READ STATS

    // get centre from stats file
    auto stats = SmallDataIO::read(m_stats + ".dat");

    // to determine time of last AH output
    // subtract 'dt=current_time / int_step' until last horizon print
    // interval
    int interval = m_params.solve_interval * m_params.print_interval;
    double current_time_reduced =
        a_current_time -
        (a_int_step % interval) * (a_current_time / a_int_step);
    pout() << "current_time_reduced = " << current_time_reduced << std::endl;
    const double time_epsilon = SmallDataIO::get_default_coords_epsilon();

    if (stats.size() == 0)
    { // case when it never ran the AH or the file doesn't exist
        pout() << "Empty stats file '" << m_stats
               << "'. Assuming AH was never found and PETSc never converged."
               << std::endl;
        m_has_been_found = false;
        m_converged = false;
    }
    else
    {
        int idx = stats[0].size();
        while (stats[0][--idx] > current_time_reduced + time_epsilon)
            ;

        if (stats[0][idx] < current_time_reduced - time_epsilon)
        { // case when the PETSc stopped converging and stopped printing
            // so there is a mismatch with the last time in the file
            pout() << "Last time step not found in '" << m_stats
                   << "'. Assuming AH was found and stopped converging."
                   << std::endl;
            m_has_been_found = true;
            m_converged = false;
        }
        else if (std::isnan(stats[1][idx]))
        { // case when the simulation was still going without ever
            // having converged
            pout()
                << "Last time step is NAN in '" << m_stats
                << "'. Assuming AH was never found and PETSc never converged."
                << std::endl;
            m_has_been_found = false;
            m_converged = false;
        }
        else
        { // case when the AH was found and there is a coordinate file
            pout() << "Last time step found in '" << m_stats
                   << "'. Reading coordinates file." << std::endl;
            m_has_been_found = true;
            m_converged = true;
        }
    }

    // if not converged, leave
    // center will remain as the initial guessed center
    // and solve() will reset coords to initial guess
    if (!m_converged)
        return;

    int rows = stats[0].size();
    int cols = stats.size();

    // look for stats line with time 'current_time_reduced', as there may be
    // further output messed up in the file after the last checkpoint was
    // made
    int i = rows;
    while (i > 0 && abs(stats[0][--i] - current_time_reduced) > time_epsilon)
        ;
    CH_assert(i > 0); // this means it was not found in the file

    pout() << "Setting center from stats file '" << m_stats << "'."
           << std::endl;

    m_center[0] = stats[cols - 3][i];
    m_center[1] = stats[cols - 2][i];
    m_center[2] = stats[cols - 1][i];
    m_center_old = m_center;
    set_center(m_center);

    // READ COORDINATES

    std::stringstream ss;
    std::string coords_filename = SmallDataIO::get_new_filename(
        m_coords, a_current_time / static_cast<double>(a_int_step),
        current_time_reduced);
    auto coords = SmallDataIO::read(coords_filename);
    int coords_ncols = coords.size();

    if (AHFinder::is_rank_active())
    {
        CH_assert(coords[0].size() == m_num_global_u * m_num_global_v);

        pout() << "Setting Initial Guess to previous file '" << coords_filename
               << "' found." << std::endl;

        // read PETSc array to 'f'
        dmda_arr_t f;
        DMDAVecGetArray(m_dmda, m_snes_soln, &f);

#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
#endif
        {
            for (int u = m_umin; u < m_umax; ++u)
            {
                int idx_global = v * m_num_global_u + u;
                int idx_local = (v - m_vmin) * m_nu + (u - m_umin);
                // m_u[idx_local] = coords[0][idx_global]; // unchanged as
                // 'write_coords_to_file' preserves order
                m_F[idx_local] = sqrt(coords[coords_ncols - 1][idx_global]);

#if CH_SPACEDIM == 3
                // m_v[idx_local] = coords[1][idx_global]; // unchanged as
                // 'write_coords_to_file' preserves order
                f[v][u] = m_F[idx_local];
#elif CH_SPACEDIM == 2
                f[u] = m_F[idx_local];
#endif
            }
        }

        // write PETSc array back
        DMDAVecRestoreArray(m_dmda, m_snes_soln, &f);
    }
}

template <typename AHGeom>
void ApparentHorizon<AHGeom>::write_coords_to_file(double a_dt, double a_time,
                                                   double a_restart_time) const
{
    CH_TIME("ApparentHorizon::write_coords_to_file");

    CH_assert(a_dt != 0); // Check if time was set!!

    SmallDataIO coords_file(m_coords, a_dt, a_time, a_restart_time,
                            SmallDataIO::NEW); // always write from scratch

    coords_file.write_header_line(m_geom.get_labels(), "");

    if (!get_converged())
    {
        coords_file.write_header_line({"Horizon NOT found."}, "");
        return;
    }

    //////////////////
    // method:
    // 1) all PESc ranks send their coordinates to rank 0
    // (non-PETSc processes will do nothing)
    // 2) rank 0 receives and stores all coordinates
    // 3) non-blocking sends finish the communication ('MPI_Wait')
    // 4) rank 0 writes all the information
    // Note: this does NOT assume that rank 0 is part of the PETSc communicator!
    // (even though for now the PETSc processes always include rank 0)
    //////////////////

    int rank = 0;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#endif

    // not needed for non-PETSc ranks, but they will be waiting anyway
    MPI_Request *requests = nullptr; // needed for non-blocking communication
    double *coords = nullptr;        // needed for non-blocking communication

    // step 1
    if (AHFinder::is_rank_active())
    {

        int local_total = (m_vmax - m_vmin) * (m_umax - m_umin);
        requests = new MPI_Request[local_total];
        coords = new double[local_total * CH_SPACEDIM];

        int idx = 0;
#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
#endif
        {
            for (int u = m_umin; u < m_umax; ++u)
            {
#if CH_SPACEDIM == 3
                coords[idx * CH_SPACEDIM] = m_u[idx];
                coords[idx * CH_SPACEDIM + 1] = m_v[idx];
                coords[idx * CH_SPACEDIM + 2] = m_F[idx] * m_F[idx];
#elif CH_SPACEDIM == 2
                coords[idx * CH_SPACEDIM] = m_v[idx];
                coords[idx * CH_SPACEDIM + 1] = m_F[idx] * m_F[idx];
#endif

#ifdef CH_MPI
                // all processes send their 'coords' to rank 0, who receives and
                // writes everything
                // (to simplify, rank 0 also sends it to itself, otherwise we
                // wouldn't know which tags rank 0 has)
                // sends are tagged by global index, so that receives can be
                // unique and writes indexed correctly
                int idx_global = v * m_num_global_u + u;
                MPI_Issend(&coords[idx * CH_SPACEDIM], CH_SPACEDIM, MPI_DOUBLE,
                           0, idx_global, Chombo_MPI::comm, &(requests[idx]));
#else
                // only rank 0 is active -> serial
                coords_file.write_data_line(
                    std::vector<double>(coords, coords + CH_SPACEDIM));
#endif

                idx++;
            }
        }
    }

    // step 2
#ifdef CH_MPI
    int total = m_num_global_u * m_num_global_v;
    double temp[total * CH_SPACEDIM];
    if (rank == 0)
    {
        for (unsigned i = 0; i < total; ++i)
        {
            MPI_Status status;
            MPI_Recv(&temp[i * CH_SPACEDIM], CH_SPACEDIM, MPI_DOUBLE,
                     MPI_ANY_SOURCE, i, Chombo_MPI::comm, &status);
        }
    }

    // step 3
    if (AHFinder::is_rank_active())
    {
        int idx = 0;
#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
#endif
        {
            for (int u = m_umin; u < m_umax; ++u)
            {
                MPI_Wait(&(requests[idx]), MPI_STATUS_IGNORE);
                idx++;
            }
        }

        delete[] requests; // free allocated memory
        delete[] coords;
    }

    // step 4
    if (rank == 0)
    {
        for (unsigned i = 0; i < total; ++i)
        {
            coords_file.write_data_line(std::vector<double>(
                &temp[i * CH_SPACEDIM], &temp[i * CH_SPACEDIM] + CH_SPACEDIM));
        }
    }
#endif
}

// ONLY FOR 3D
template <typename AHGeom>
double ApparentHorizon<AHGeom>::calculate_spin(double a_area)
{
    CH_assert(CH_SPACEDIM == 3);
    CH_TIME("ApparentHorizon::calculate_spin");

    if (!get_converged())
        return NAN;

    double equator_length;
    double tempC = 0.; // temporary, but defined outside to be in scope for
                       // non-PETSc processes

    // PETSc processes go inside 'if', others "wait" until 'if' gets to
    // 'm_geom.break_interpolation_loop()'
    if (m_geom.keep_interpolating_if_inactive())
    {

        int idx = 0;

        Vec localF;

        DMGetLocalVector(m_dmda, &localF);
        DMGlobalToLocalBegin(m_dmda, m_snes_soln, INSERT_VALUES, localF);
        DMGlobalToLocalEnd(m_dmda, m_snes_soln, INSERT_VALUES, localF);

        dmda_arr_t in;
        DMDAVecGetArray(m_dmda, localF, &in);

#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
#endif
        {
            for (int u = m_umin; u < m_umax; ++u)
            {
#if CH_SPACEDIM == 3
                m_F[(v - m_vmin) * m_nu + (u - m_umin)] = in[v][u];
#elif CH_SPACEDIM == 2
                m_F[u - m_umin] = in[u];
#endif
            }
        }

        m_geom.set_data(m_u, m_v, m_F);

        double norm2 = 0;
        double dxdv[3];
        double dxdu[3];

        int u_equator = floor(M_PI / 2.0 / m_du);

        for (int v = m_vmin; v < m_vmax; ++v)
        {
            for (int u = m_umin; u < m_umax; ++u)
            {
                Deriv deriv = diff(in, u, v);

                if (u_equator == u)
                {

                    const AHGeometryData data = m_geom.get_data(idx);

                    dxdv[0] = -pow(m_F[idx], 2) * sin(m_u[idx]) * sin(m_v[idx]);
                    +2 * m_F[idx] * deriv.dvF *sin(m_u[idx]) * cos(m_v[idx]);

                    dxdv[1] = pow(m_F[idx], 2) * sin(m_u[idx]) * cos(m_v[idx]);
                    +2 * m_F[idx] * deriv.dvF *sin(m_u[idx]) * sin(m_v[idx]);

                    dxdv[2] = 0 + 2 * m_F[idx] * deriv.dvF * cos(m_u[idx]);

                    norm2 = 0;
                    FOR2(i, j) norm2 += data.g[i][j] * dxdv[i] * dxdv[j];

                    CH_assert(norm2 >= 0.);

                    tempC += sqrt(norm2) * m_dv;
                }

                idx++;
            }
        }

        DMDAVecRestoreArray(m_dmda, localF, &in);
        DMRestoreLocalVector(m_dmda, &localF);

        m_geom.break_interpolation_loop();
    }

    // reduction across all Chombo processes (note that 'tempC' is 0 for
    // non-PETSc processes) because SmallDataIO uses rank 0 to write this
    // ensures rank 0 will have 'tempC' (even though for now the PETSc processes
    // always include rank 0)
#ifdef CH_MPI
    MPI_Allreduce(&tempC, &equator_length, 1, MPI_DOUBLE, MPI_SUM,
                  Chombo_MPI::comm);
#else // serial
    equator_length = tempC;
#endif

    double factor =
        ((2. * M_PI * a_area / (equator_length * equator_length)) - 1.);
    double spin =
        (factor > 1. ? 0
                     : sqrt(1. - factor * factor)); // factor>1 means numerical
                                                    // error with spin as 0

    pout() << "spin = " << spin << endl;

    return spin;
}

template <typename AHGeom>
double ApparentHorizon<AHGeom>::calculate_coord_area()
{
    CH_assert(CH_SPACEDIM == 3);
    CH_TIME("ApparentHorizon::calculate_coord_area");

    if (!get_converged())
        return NAN;

    double coord_area;
    double temp =
        0.; // temporary, but defined outside to exist for non-PETSc processes
    const std::string tag = "ApparentHorizon::calculate_coord_area: ";
    // PETSc processes go inside 'if', others "wait" until 'if' gets to
    // 'm_geom.break_interpolation_loop()'
    if (m_geom.keep_interpolating_if_inactive())
    {

        int idx = 0;

        Vec localF;

        DMGetLocalVector(m_dmda, &localF);
        DMGlobalToLocalBegin(m_dmda, m_snes_soln, INSERT_VALUES, localF);
        DMGlobalToLocalEnd(m_dmda, m_snes_soln, INSERT_VALUES, localF);

        dmda_arr_t in;
        DMDAVecGetArray(m_dmda, localF, &in);

#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
#endif
        {
            for (int u = m_umin; u < m_umax; ++u)
            {
#if CH_SPACEDIM == 3
                m_F[(v - m_vmin) * m_nu + (u - m_umin)] = in[v][u];
#elif CH_SPACEDIM == 2
                m_F[u - m_umin] = in[u];
#endif
            }
        }

        m_geom.set_data(m_u, m_v, m_F);

#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
#endif
        {
            for (int u = m_umin; u < m_umax; ++u)
            {

                const AHGeometryData data = m_geom.get_data(idx);

                // Calculate Jacobian matrix for transformation from Cartesian
                // to (f,u,v) coords
                Tensor<2, double> Jac;
                FOR1(k)
                {
                    Jac[0][k] = data.dxdf[k];
                    Jac[1][k] = data.dxdu[k];
                    Jac[2][k] = data.dxdv[k];
                }
                // auto det = TensorAlgebra::compute_determinant(Jac);
                Deriv deriv = diff(in, u, v);

                // Now do the coordinate transformation
                Tensor<2, double> g_spherical = {0.};
                Tensor<2, double> g_flat = {0.};
                g_flat[0][0] = g_flat[1][1] = g_flat[2][2] = 1.0;
                FOR4(i, j, k, l)
                {
                    g_spherical[i][j] += Jac[i][k] * Jac[j][l] * g_flat[k][l];
                }

                // Construct the 2-metric on the horizon in (u,v) coords
                // i.e. substitute df = (df/du)du + (df/dv)dv
                // into the spherical metric
                Tensor<2, double, 2> g_horizon = {0.};
                g_horizon[0][0] = g_spherical[1][1] +
                                  g_spherical[0][0] * deriv.duF * deriv.duF +
                                  2.0 * g_spherical[0][1] * deriv.duF;
                g_horizon[1][1] = g_spherical[2][2] +
                                  g_spherical[0][0] * deriv.dvF * deriv.dvF +
                                  2.0 * g_spherical[0][2] * deriv.dvF;
                g_horizon[0][1] = g_horizon[1][0] =
                    g_spherical[1][2] +
                    g_spherical[0][0] * deriv.duF * deriv.dvF +
                    g_spherical[0][1] * deriv.dvF +
                    g_spherical[0][2] * deriv.duF;

                double det = TensorAlgebra::compute_determinant(g_horizon);

                if (!m_periodic_u && (u == 0 || u == m_num_global_u - 1))
                {
                }
#if CH_SPACEDIM == 3
                else if (!m_periodic_v && (v == 0 || v == m_num_global_v - 1))
                {
                    temp += 0.5 * sqrt(det) * m_du * m_dv;
                }
#endif
                else
                {
#if CH_SPACEDIM == 3
                    temp += sqrt(det) * m_du * m_dv;
#endif
                }
                idx++;
            }
        }

        DMDAVecRestoreArray(m_dmda, localF, &in);
        DMRestoreLocalVector(m_dmda, &localF);

        m_geom.break_interpolation_loop();
    }

// reduction across all Chombo processes (note that 'area' is 0 for non-PETSc
// processes) because SmallDataIO uses rank 0 to write this ensures rank 0 will
// have the area (even though for now the PETSc processes always include rank 0)
#ifdef CH_MPI
    // we want all the processes to have the area so that all use it in the
    // spin calculation (which in reality is needed not for the output files,
    // but only for the pout() prints)
    MPI_Allreduce(&temp, &coord_area, 1, MPI_DOUBLE, MPI_SUM, Chombo_MPI::comm);
#else // serial
    coord_area = temp;
#endif

    pout() << "coord area = " << coord_area << endl;

    return coord_area;
}

template <typename AHGeom> double ApparentHorizon<AHGeom>::calculate_area()
{
    CH_assert(CH_SPACEDIM == 3);
    CH_TIME("ApparentHorizon::calculate_area");

    if (!get_converged())
        return NAN;

    double area;
    double temp =
        0.; // temporary, but defined outside to exist for non-PETSc processes
    const std::string tag = "ApparentHorizon::calculate_area: ";
    // PETSc processes go inside 'if', others "wait" until 'if' gets to
    // 'm_geom.break_interpolation_loop()'
    if (m_geom.keep_interpolating_if_inactive())
    {

        int idx = 0;

        Vec localF;

        DMGetLocalVector(m_dmda, &localF);
        DMGlobalToLocalBegin(m_dmda, m_snes_soln, INSERT_VALUES, localF);
        DMGlobalToLocalEnd(m_dmda, m_snes_soln, INSERT_VALUES, localF);

        dmda_arr_t in;
        DMDAVecGetArray(m_dmda, localF, &in);

#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
#endif
        {
            for (int u = m_umin; u < m_umax; ++u)
            {
#if CH_SPACEDIM == 3
                m_F[(v - m_vmin) * m_nu + (u - m_umin)] = in[v][u];
#elif CH_SPACEDIM == 2
                m_F[u - m_umin] = in[u];
#endif
            }
        }

        m_geom.set_data(m_u, m_v, m_F);

#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
#endif
        {
            for (int u = m_umin; u < m_umax; ++u)
            {

                double area_element = 0;
                const AHGeometryData data = m_geom.get_data(idx);

                // Calculate Jacobian matrix for transformation from Cartesian
                // to (f,u,v) coords
                Tensor<2, double> Jac;
                FOR1(k)
                {
                    Jac[0][k] = data.dxdf[k];
                    Jac[1][k] = data.dxdu[k];
                    Jac[2][k] = data.dxdv[k];
                }
                // auto det = TensorAlgebra::compute_determinant(Jac);
                Deriv deriv = diff(in, u, v);

                // Now do the coordinate transformation
                Tensor<2, double> g_spherical = {0.};
                FOR4(i, j, k, l)
                {
                    g_spherical[i][j] += Jac[i][k] * Jac[j][l] * data.g[k][l];
                }

                // Construct the 2-metric on the horizon in (u,v) coords
                // i.e. substitute df = (df/du)du + (df/dv)dv
                // into the spherical metric
                Tensor<2, double, 2> g_horizon = {0.};
                g_horizon[0][0] = g_spherical[1][1] +
                                  g_spherical[0][0] * deriv.duF * deriv.duF +
                                  2.0 * g_spherical[0][1] * deriv.duF;
                g_horizon[1][1] = g_spherical[2][2] +
                                  g_spherical[0][0] * deriv.dvF * deriv.dvF +
                                  2.0 * g_spherical[0][2] * deriv.dvF;
                g_horizon[0][1] = g_horizon[1][0] =
                    g_spherical[1][2] +
                    g_spherical[0][0] * deriv.duF * deriv.dvF +
                    g_spherical[0][1] * deriv.dvF +
                    g_spherical[0][2] * deriv.duF;

                double det = TensorAlgebra::compute_determinant(g_horizon);

                if (!m_periodic_u && (u == 0 || u == m_num_global_u - 1))
                {
                }
#if CH_SPACEDIM == 3
                else if (!m_periodic_v && (v == 0 || v == m_num_global_v - 1))
                {
                    temp += 0.5 * sqrt(det) * m_du * m_dv;
                }
#endif
                else
                {
#if CH_SPACEDIM == 3
                    /*
                    #include "INTFunction.inc"
                                        double rr = pow(m_F[idx], 2);
                                        double drdu = 2 * m_F[idx] * deriv.duF;
                                        double drdv = 2 * m_F[idx] * deriv.dvF;
                                        double coordinate_element =
                                            rr * rr *
                                            (rr * rr * pow(sin(m_u[idx]), 2) +
                    drdv * drdv + drdu * drdu * pow(sin(m_u[idx]), 2));

                                        pout() << tag << "coordinate_element("
                    << m_u[idx] << ", "
                                               << m_v[idx] << ") = " <<
                    coordinate_element << "\n"; pout() << tag << "area_alement("
                    << m_u[idx] << ", "
                                              << m_v[idx] << ") = " <<
                    area_element << std::endl;

                                        temp += sqrt(area_element) *
                    sqrt(coordinate_element) * m_du * m_dv;
                                        */
                    temp += sqrt(det) * m_du * m_dv;
#endif
                }
                idx++;
            }
        }

        DMDAVecRestoreArray(m_dmda, localF, &in);
        DMRestoreLocalVector(m_dmda, &localF);

        m_geom.break_interpolation_loop();
    }

// reduction across all Chombo processes (note that 'area' is 0 for non-PETSc
// processes) because SmallDataIO uses rank 0 to write this ensures rank 0 will
// have the area (even though for now the PETSc processes always include rank 0)
#ifdef CH_MPI
    // we want all the processes to have the area so that all use it in the
    // spin calculation (which in reality is needed not for the output files,
    // but only for the pout() prints)
    MPI_Allreduce(&temp, &area, 1, MPI_DOUBLE, MPI_SUM, Chombo_MPI::comm);
#else // serial
    area = temp;
#endif

    pout() << "area = " << area << endl;

    return area;
}

template <typename AHGeom>
std::array<double, CH_SPACEDIM> ApparentHorizon<AHGeom>::calculate_center()
{
    CH_assert(CH_SPACEDIM == 3);
    CH_TIME("ApparentHorizon::calculate_center");

    if (!get_converged())
        return {NAN, NAN, NAN};

    // Method:
    // Calculate centroid of all the points, by summing the coordinates {x,y,z}
    // of all the points, which are spread across processors, and then reducing
    // them all with MPI. Finally, divide by total number of points summed, to
    // get the centroid

    std::array<double, CH_SPACEDIM> temp = {0., 0., 0.};

    int idx = 0;
    if (AHFinder::is_rank_active()) // in principle this wouldn't be needed, as
                                    // non-PETSc processes wouldn't even enter
                                    // the loop anyways
    {
#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
#endif
        {
            for (int u = m_umin; u < m_umax; ++u)
            {
                // actually transform() returns a 'std::array<double, 3>',
                // without the 'CH_SPACEDIM'
                std::array<double, CH_SPACEDIM> point =
                    m_geom.transform(m_u[idx], m_v[idx], m_F[idx]);

                for (unsigned i = 0; i < CH_SPACEDIM; ++i)
                    temp[i] += point[i];

                idx++;
            }
        }
    }

    std::array<double, CH_SPACEDIM> center;
    int idx_sum;

#ifdef CH_MPI
    MPI_Allreduce(&temp, &center, CH_SPACEDIM, MPI_DOUBLE, MPI_SUM,
                  Chombo_MPI::comm);
    MPI_Allreduce(&idx, &idx_sum, 1, MPI_INT, MPI_SUM, Chombo_MPI::comm);
#else // serial
    center = temp;
    idx_sum = idx;
#endif

    CH_assert(idx_sum != 0);

    for (unsigned i = 0; i < CH_SPACEDIM; ++i)
        center[i] /= idx_sum;

    set_center(center);

    pout() << "center: (" << center[0] << ", " << center[1] << ", " << center[2]
           << ")" << std::endl;

    return center;
}

template <typename AHGeom> void ApparentHorizon<AHGeom>::calculate_minmax_F()
{
    double local_max = 0;
    double local_min = 0;
    if (AHFinder::is_rank_active())
    {
        auto local_minmax = std::minmax_element(m_F.begin(), m_F.end());
        local_min = *(local_minmax.first);
        local_max = *(local_minmax.second);
    }
    double global_max = local_max;
    double global_min = local_min;

#ifdef CH_MPI
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN,
                  Chombo_MPI::comm);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX,
                  Chombo_MPI::comm);
#endif

    m_max_F = global_max;
    m_min_F = global_min;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////// PETSc stuff below ///////////////////////////
/////////////////////////////////////////////////////////////////////////

template <typename AHGeom> void ApparentHorizon<AHGeom>::initialise_PETSc()
{
    CH_TIME("ApparentHorizon::initialise_PETSc");

    if (!AHFinder::is_rank_active())
        return;

    const double lo_u = m_geom.get_domain_low()[0];
    const double hi_u = m_geom.get_domain_high()[0];
#if CH_SPACEDIM == 3 // lo_v, hi_v not used otherwise
    const double lo_v = m_geom.get_domain_low()[1];
    const double hi_v = m_geom.get_domain_high()[1];
#endif

#if PETSC_VERSION_LT(3, 5, 0)
#define DM_BOUNDARY_PERIODIC DMDA_BOUNDARY_PERIODIC
#define DM_BOUNDARY_GHOSTED DMDA_BOUNDARY_PERIODIC
#endif

#if CH_SPACEDIM == 3
    DMDACreate2d(PETSC_COMM_WORLD,
                 m_periodic_u ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                 m_periodic_v ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                 DMDA_STENCIL_BOX, m_num_global_u,
                 m_num_global_v, /* grid size (negative means that we can
                                    override it using command line) */
                 PETSC_DECIDE, PETSC_DECIDE, /* distribution */
                 1,                          /* number of degrees of freedom */
                 3, /* stencil width (each side from central point) */
                 NULL, NULL, &m_dmda);
#elif CH_SPACEDIM == 2
    DMDACreate1d(PETSC_COMM_WORLD,
                 m_periodic_u ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                 m_num_global_u, /* grid size (negative means that we can
                                    override it using command line) */
                 1,              /* number of degrees of freedom */
                 3, /* stencil width (each side from central point) */
                 NULL, &m_dmda);
#endif

    DMSetUp(m_dmda);

    // reload -> is it the same value?
    DMDAGetInfo(m_dmda, NULL, &m_num_global_u, &m_num_global_v, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    m_du =
        (hi_u - lo_u) / (m_periodic_u ? m_num_global_u : (m_num_global_u - 1));

#if CH_SPACEDIM == 3
    m_dv =
        (hi_v - lo_v) / (m_periodic_v ? m_num_global_v : (m_num_global_v - 1));
#endif

    DMDAGetCorners(m_dmda, &m_umin, &m_vmin, NULL, &m_nu, &m_nv, NULL);
    m_umax = m_umin + m_nu;
    m_vmax = m_vmin + m_nv;

#if CH_SPACEDIM == 3
    int vec_size = m_nu * m_nv;
#elif CH_SPACEDIM == 2
    int vec_size = m_nu;
#endif

    m_u.resize(vec_size);
    m_v.resize(vec_size);
    m_F.resize(vec_size);

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
#endif
    {
        for (int u = m_umin; u < m_umax; ++u)
        {
#if CH_SPACEDIM == 3
            m_u[(v - m_vmin) * m_nu + (u - m_umin)] = lo_u + u * m_du;
            m_v[(v - m_vmin) * m_nu + (u - m_umin)] = lo_v + v * m_dv;
#elif CH_SPACEDIM == 2
            m_u[u - m_umin] = lo_u + u * m_du;
            m_v[u - m_umin] = 0;
#endif
        }
    }

    DMCreateGlobalVector(m_dmda, &m_snes_soln);
    PetscObjectSetName((PetscObject)m_snes_soln, "F");

    VecDuplicate(m_snes_soln, &m_snes_rhs);
    PetscObjectSetName((PetscObject)m_snes_rhs, "expansion");

    VecDuplicate(m_snes_soln, &m_snes_guu);
    PetscObjectSetName((PetscObject)m_snes_guu, "guu");

    VecDuplicate(m_snes_soln, &m_snes_gvv);
    PetscObjectSetName((PetscObject)m_snes_gvv, "gvv");

    VecDuplicate(m_snes_soln, &m_snes_guv);
    PetscObjectSetName((PetscObject)m_snes_guv, "guv");

    VecDuplicate(m_snes_soln, &m_snes_gww);
    PetscObjectSetName((PetscObject)m_snes_gww, "gww");

    VecDuplicate(m_snes_soln, &m_snes_lapse);
    PetscObjectSetName((PetscObject)m_snes_lapse, "lapse");

#if PETSC_VERSION_GE(3, 5, 0)
    DMSetMatType(m_dmda, MATAIJ);
    DMCreateMatrix(m_dmda, &m_snes_jac);
    MatSetOption(m_snes_jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#else
    DMCreateMatrix(m_dmda, MATAIJ, &m_snes_jac);
    MatSetOption(m_snes_jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#endif

    SNESCreate(PETSC_COMM_WORLD, &m_snes);
    SNESSetFunction(m_snes, m_snes_rhs, &Petsc_form_function, this);
    SNESSetJacobian(m_snes, m_snes_jac, m_snes_jac, &Petsc_form_jacobian, this);
    SNESMonitorSet(m_snes, &Petsc_SNES_monitor, this, NULL);

    SNESSetFromOptions(m_snes);

    m_geom.set_SNES(m_snes);
    m_geom_plus.set_SNES(m_snes);
    m_geom_minus.set_SNES(m_snes);

    SNESType snes_type;
    SNESGetType(m_snes, &snes_type);
    double snes_atol, snes_rtol, snes_stol;
    int snes_maxit, snes_maxf;
    SNESGetTolerances(m_snes, &snes_atol, &snes_rtol, &snes_stol, &snes_maxit,
                      &snes_maxf);
    KSP snes_ksp;
    SNESGetKSP(m_snes, &snes_ksp);
    KSPType ksp_type;
    KSPGetType(snes_ksp, &ksp_type);
    double ksp_rtol, ksp_abstol, ksp_dtol;
    int ksp_maxits;
    KSPGetTolerances(snes_ksp, &ksp_rtol, &ksp_abstol, &ksp_dtol, &ksp_maxits);

    pout() << "-------------------------------------\n";
    pout() << "SNES Options:\n";
    pout() << "Type: " << snes_type << "\n";
    pout() << "atol = " << snes_atol << ", rtol = " << snes_rtol
           << ", stol = " << snes_stol << ",\n";
    pout() << "maxit = " << snes_maxit << ", maxf = " << snes_maxf << "\n";
    pout() << "-------------------------------------\n";
    pout() << "KSP Options:\n";
    pout() << "Type: " << ksp_type << "\n";
    pout() << "rtol = " << ksp_rtol << ", abstol = " << ksp_abstol
           << ", dtol = " << ksp_dtol << ", maxits = " << ksp_maxits << "\n";
    pout() << "-------------------------------------" << std::endl;
}

template <typename AHGeom> void ApparentHorizon<AHGeom>::finalise_PETSc()
{
    if (!AHFinder::is_rank_active())
        return;

    SNESDestroy(&m_snes);
    VecDestroy(&m_snes_soln);
    VecDestroy(&m_snes_rhs);
    MatDestroy(&m_snes_jac);
    VecDestroy(&m_snes_guu);
    VecDestroy(&m_snes_gvv);
    VecDestroy(&m_snes_guv);
    VecDestroy(&m_snes_gww);
    VecDestroy(&m_snes_lapse);
    DMDestroy(&m_dmda);
}

template <typename AHGeom>
typename ApparentHorizon<AHGeom>::Deriv
ApparentHorizon<AHGeom>::diff(dmda_arr_t in, int u, int v)
{
    CH_TIME("ApparentHorizon::diff");

    Deriv out;
#include "AHStencil.inc"

#if CH_SPACEDIM == 3

    // d/du
    for (int j = 0; j < DWIDTH; ++j)
    {
        out.duF +=
            out.du_weights[j] * in[v][u + out.du_stencil_start + j] / m_du;

        // d2/dudv
        for (int k = 0; k < DWIDTH; ++k)
        {
            out.dudvF +=
                out.du_weights[j] * out.dv_weights[k] *
                in[v + out.dv_stencil_start + k][u + out.du_stencil_start + j] /
                (m_du * m_dv);
        }
    }

    // d/dv
    for (int j = 0; j < DWIDTH; ++j)
    {
        out.dvF +=
            out.dv_weights[j] * in[v + out.dv_stencil_start + j][u] / m_dv;
    }

    // d2/du2
    for (int j = 0; j < DDWIDTH; ++j)
    {
        out.duduF += out.dudu_weights[j] *
                     in[v][u + out.dudu_stencil_start + j] / (m_du * m_du);
    }

    // d2/dv2
    for (int j = 0; j < DDWIDTH; ++j)
    {
        out.dvdvF += out.dvdv_weights[j] *
                     in[v + out.dvdv_stencil_start + j][u] / (m_dv * m_dv);
    }

#elif CH_SPACEDIM == 2

    // d/du
    for (int j = 0; j < DWIDTH; ++j)
    {
        out.duF += out.du_weights[j] * in[u + out.du_stencil_start + j] / m_du;
    }

    // d2/du2
    for (int j = 0; j < DDWIDTH; ++j)
    {
        out.duduF += out.dudu_weights[j] * in[u + out.dudu_stencil_start + j] /
                     (m_du * m_du);
    }

#endif

    return out;
}

template <typename AHGeom>
void ApparentHorizon<AHGeom>::form_function(Vec F, Vec Rhs)
{
    CH_TIME("ApparentHorizon::form_function");

    // Scatter ghost cells
    Vec localF;
    DMGetLocalVector(m_dmda, &localF);
    DMGlobalToLocalBegin(m_dmda, F, INSERT_VALUES, localF);
    DMGlobalToLocalEnd(m_dmda, F, INSERT_VALUES, localF);

    dmda_arr_t in;
    DMDAVecGetArray(m_dmda, localF, &in);

    dmda_arr_t out;
    DMDAVecGetArray(m_dmda, Rhs, &out);

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
#endif
    {
        for (int u = m_umin; u < m_umax; ++u)
        {
#if CH_SPACEDIM == 3
            m_F[(v - m_vmin) * m_nu + (u - m_umin)] = in[v][u] + eps;
#elif CH_SPACEDIM == 2
            m_F[u - m_umin] = in[u] + eps;
#endif
        }
    }

    m_geom_plus.set_data(m_u, m_v, m_F);

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
#endif
    {
        for (int u = m_umin; u < m_umax; ++u)
        {
#if CH_SPACEDIM == 3
            m_F[(v - m_vmin) * m_nu + (u - m_umin)] = in[v][u] - eps;
#elif CH_SPACEDIM == 2
            m_F[u - m_umin] = in[u] - eps;
#endif
        }
    }

    m_geom_minus.set_data(m_u, m_v, m_F);

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
#endif
    {
        for (int u = m_umin; u < m_umax; ++u)
        {
#if CH_SPACEDIM == 3
            m_F[(v - m_vmin) * m_nu + (u - m_umin)] = in[v][u];
#elif CH_SPACEDIM == 2
            m_F[u - m_umin] = in[u];
#endif
        }
    }

    m_geom.set_data(m_u, m_v, m_F);

    int idx = 0;

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
    {
#elif CH_SPACEDIM == 2
    {
        int v = 0;
#endif

        for (int u = m_umin; u < m_umax; ++u)
        {

#if CH_SPACEDIM == 3
            double &_out = out[v][u];
#elif CH_SPACEDIM == 2
            double &_out = out[u];
#endif

            Deriv deriv = diff(in, u, v);

            if (!m_periodic_u && (u == 0 || u == m_num_global_u - 1))
            {
                _out = deriv.duF;
            }

#if CH_SPACEDIM == 3
            else if (!m_periodic_v && (v == 0 || v == m_num_global_v - 1))
            {
                _out = deriv.dvF;
            }
#endif

            else
            {
                const AHGeometryData data = m_geom.get_data(idx);

                double expansion = 0;
#include "AHFunction.inc"

                _out = expansion;
            }

            ++idx;
        }
    }

    DMDAVecRestoreArray(m_dmda, localF, &in);
    DMDAVecRestoreArray(m_dmda, Rhs, &out);
    DMRestoreLocalVector(m_dmda, &localF);
}

template <typename AHGeom>
void ApparentHorizon<AHGeom>::form_jacobian(Vec F, Mat J)
{
    CH_TIME("ApparentHorizon::form_jacobian");

    // Scatter ghost cells
    Vec localF;
    DMGetLocalVector(m_dmda, &localF);
    DMGlobalToLocalBegin(m_dmda, F, INSERT_VALUES, localF);
    DMGlobalToLocalEnd(m_dmda, F, INSERT_VALUES, localF);

    dmda_arr_t in;
    DMDAVecGetArray(m_dmda, localF, &in);

    int idx = 0;

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
    {
#elif CH_SPACEDIM == 2
    {
        int v = 0;
#endif

        for (int u = m_umin; u < m_umax; ++u)
        {
            MatStencil row[1] = {0};
            row[0].i = u;
            row[0].j = v;

            if (!m_periodic_u && (u == 0 || u == m_num_global_u - 1))
            {
                MatStencil col[DWIDTH] = {0};
                double val[DWIDTH] = {0};

                const Deriv deriv = diff(in, u, v);

                for (int a = 0; a < DWIDTH; ++a)
                {
                    col[a].i = u + deriv.du_stencil_start + a;
                    col[a].j = v;
                    val[a] = deriv.du_weights[a];
                }

                MatSetValuesStencil(J, 1, row, DWIDTH, col, val, INSERT_VALUES);
            }

#if CH_SPACEDIM == 3
            else if (!m_periodic_v && (v == 0 || v == m_num_global_v - 1))
            {
                MatStencil col[DWIDTH] = {0};
                double val[DWIDTH] = {0};

                const Deriv deriv = diff(in, u, v);

                for (int b = 0; b < DWIDTH; ++b)
                {
                    col[b].i = u;
                    col[b].j = v + deriv.dv_stencil_start + b;
                    val[b] = deriv.dv_weights[b];
                }

                MatSetValuesStencil(J, 1, row, DWIDTH, col, val, INSERT_VALUES);
            }
#endif

            else
            {

#if CH_SPACEDIM == 3
                const int NVAL = DDWIDTH * DDWIDTH;
#elif CH_SPACEDIM == 2
                const int NVAL = DDWIDTH;
#endif

                MatStencil col[NVAL] = {0};
                double val[NVAL] = {0};

                const Deriv deriv_default = diff(in, u, v);

#if CH_SPACEDIM == 3
                for (int b = 0; b < DDWIDTH; ++b)
                {
#elif CH_SPACEDIM == 2
                {
                    int b = 0;
#endif

                    for (int a = 0; a < DDWIDTH; ++a)
                    {
                        col[b * DDWIDTH + a].i =
                            u + deriv_default.dudu_stencil_start + a;
                        col[b * DDWIDTH + a].j =
                            v + deriv_default.dvdv_stencil_start + b;
                    }
                }

                // "local" jacobian term
                {
#if CH_SPACEDIM == 3
                    double &_in = in[v][u];
#elif CH_SPACEDIM == 2
                    double &_in = in[u];
#endif

                    // "plus" perturbation
                    double expansionPlus = 0.0;
                    {
                        double in_old = _in;
                        _in += eps;

                        const Deriv deriv = diff(in, u, v);
                        const AHGeometryData data = m_geom_plus.get_data(idx);

                        double expansion = 0.0;
#include "AHFunction.inc"

                        expansionPlus = expansion;
                        _in = in_old;
                    }

                    // "minus" perturbation
                    double expansionMinus = 0.0;
                    {
                        double in_old = _in;
                        _in -= eps;

                        const Deriv deriv = diff(in, u, v);
                        const AHGeometryData data = m_geom_minus.get_data(idx);

                        double expansion = 0.0;
#include "AHFunction.inc"

                        expansionMinus = expansion;
                        _in = in_old;
                    }

                    const int a = -deriv_default.dudu_stencil_start;

#if CH_SPACEDIM == 3
                    const int b = -deriv_default.dvdv_stencil_start;
#elif CH_SPACEDIM == 2
                    const int b = 0;
#endif

                    val[b * DDWIDTH + a] =
                        (expansionPlus - expansionMinus) / (2 * eps);
                }

                // "stencil" jacobian terms
                {
#if CH_SPACEDIM == 3
                    for (int b = 0; b < DDWIDTH; ++b)
                    {
#elif CH_SPACEDIM == 2
                    {
                        int b = 0;
#endif

                        for (int a = 0; a < DDWIDTH; ++a)
                        {
                            if (a == -deriv_default.dudu_stencil_start &&
                                b == -deriv_default.dvdv_stencil_start)
                            {
                                continue;
                            }

                            const int uu = col[b * DDWIDTH + a].i;
                            const int vv = col[b * DDWIDTH + a].j;

#if CH_SPACEDIM == 3
                            double &_in = in[vv][uu];
#elif CH_SPACEDIM == 2
                            double &_in = in[uu];
#endif

                            // "plus" perturbation
                            double expansionPlus = 0.0;
                            {
                                double in_old = _in;
                                _in += eps;

                                const Deriv deriv = diff(in, u, v);
                                const AHGeometryData data =
                                    m_geom.get_data(idx);

                                double expansion = 0.0;
#include "AHFunction.inc"

                                expansionPlus = expansion;
                                _in = in_old;
                            }

                            // "minus" perturbation
                            double expansionMinus = 0.0;
                            {
                                double in_old = _in;
                                _in -= eps;

                                const Deriv deriv = diff(in, u, v);
                                const AHGeometryData data =
                                    m_geom.get_data(idx);

                                double expansion = 0.0;
#include "AHFunction.inc"

                                expansionMinus = expansion;
                                _in = in_old;
                            }

                            val[b * DDWIDTH + a] =
                                (expansionPlus - expansionMinus) / (2 * eps);
                        }
                    }
                }

                MatSetValuesStencil(J, 1, row, NVAL, col, val, INSERT_VALUES);
            }

            ++idx;
        }
    }

    DMDAVecRestoreArray(m_dmda, localF, &in);
    DMRestoreLocalVector(m_dmda, &localF);

    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
}

#endif /* _APPARENTHORIZON_IMPL_HPP_ */
