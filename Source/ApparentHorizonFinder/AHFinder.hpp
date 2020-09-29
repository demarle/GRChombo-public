/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef USE_AHFINDER

#ifndef _AHFINDER_HPP_
#define _AHFINDER_HPP_

#ifdef CH_MPI
#include <mpi.h>
#endif

#include "SPMD.H" //Chombo_MPI::comm
#include <petsc.h>

//! Class to control PETSc MPI sub-communicator

class BHAMR;

namespace AHFinder
{

struct params
{
    int num_ranks; //!< number of ranks for PETSc sub-communicator
    int num_points_u, num_points_v; //!< number of points for 2D coordinate grid
    int solve_interval;             //!< same as checkpoint_interval, for
                                    //!< ApparentHorizon::solve
    int print_interval;             //!< same as solve_interval, but for prints
    int max_failed_convergences;    //!< the number of failed convergences
                                    //!< before solving is stopped (0 = inf)
    bool updateCenter; //!< whether or not to update the center during the
                       //!< evolution (set to false if you know it won't move)
    int level;         //!< the level to do apparent horizon finding on
};

void add_ah(
    BHAMR &a_bhamr,
    const std::array<double, CH_SPACEDIM>
        &a_center,          //!< Initial guess for center
    double a_initial_guess, //!< Initial guess for radius (or whatever
                            //!< coordinate you're solving for)
    const params &a_params, //!< set of AH parameters
    const std::string &a_stats = "stats",  //!< name for stats file with
                                           //!< area, spin and AH center
    const std::string &a_coords = "coords" //!< name for coords file with AH
                                           //!< coordinates at each time step)
);

void add_ah_merger(
    BHAMR &a_bhamr, int ah1, int ah2,
    const params &a_params,                //!< set of AH parameters
    const std::string &a_stats = "stats",  //!< name for stats file with
                                           //!< area, spin and AH center
    const std::string &a_coords = "coords" //!< name for coords file with AH
                                           //!< coordinates at each time step)
);

//! Find AH; Calculate area and spin; Update center; Print outputs
void solve(BHAMR &a_bhamr, double a_dt, double a_time, double a_restart_time);

/////////////////////////////////////////////////////////
// PETSc control methods
/////////////////////////////////////////////////////////

extern bool m_initialized; //!< is initialized?

#ifdef CH_MPI
extern MPI_Group m_mpi_group; //!< set of MPI ranks for PETSc
extern MPI_Comm m_mpi_comm;   //!< MPI sub-communicator
#endif

//! define number of ranks of PETSc sub-communicator
void set_num_ranks(int a_num_ranks);

//! initialize PETSc and its MPI sub-communicator
PetscErrorCode initialize(int a_num_ranks, int argc = 0,
                          char *argv[] = nullptr);

//! finalize PETSc
PetscErrorCode finalize();

//! true if part of PETSc MPI sub-communicator
bool is_rank_active();

}; // namespace AHFinder

#endif /* _AHFINDER_HPP_ */

#endif /* USE_AHFINDER */
