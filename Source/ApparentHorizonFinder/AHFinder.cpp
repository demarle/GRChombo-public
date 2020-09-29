/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef USE_AHFINDER

#include "AHFinder.hpp"
#include "BHAMR.hpp"

void AHFinder::add_ah(BHAMR &a_bhamr,
                      const std::array<double, CH_SPACEDIM> &a_center,
                      double a_initial_guess, const AHFinder::params &a_params,
                      const std::string &a_stats, const std::string &a_coords)
{
    if (!AHFinder::m_initialized)
        AHFinder::initialize(a_params.num_ranks);

    // determine how many AH there are already
    int num_ah = a_bhamr.m_apparent_horizons.size();
    std::stringstream ss;
    ss << "_AH" << (num_ah + 1);

    AHInterpolation<AHSphericalCoords> geom(AHSphericalCoords(),
                                            a_bhamr.m_interpolator);

    a_bhamr.m_apparent_horizons.push_back(
        new ApparentHorizon<AHInterpolation<AHSphericalCoords>>(
            geom, a_center, a_initial_guess, a_params, a_stats + ss.str(),
            a_coords + ss.str() + "_"));
}

void AHFinder::add_ah_merger(
    BHAMR &a_bhamr, int ah1, int ah2,
    const AHFinder::params &a_params, //!< set of AH parameters
    const std::string &a_stats,       //!< name for stats file with
                                      //!< area, spin and AH center
    const std::string &a_coords //!< name for coords file with AH coordinates at
                                //!< each time step)
)
{
    int num_ah = a_bhamr.m_apparent_horizons.size();
    CH_assert(ah1 >= 0 && ah1 < num_ah);
    CH_assert(ah2 >= 0 && ah2 < num_ah);
    CH_assert(ah1 != ah2);

    auto AH1 = a_bhamr.m_apparent_horizons[ah1];
    auto AH2 = a_bhamr.m_apparent_horizons[ah2];

    // set as initial guess for merged AH the sum of the other two
    double initial_guess_merged =
        AH1->get_initial_guess() + AH2->get_initial_guess();

    std::array<double, CH_SPACEDIM> center1 = AH1->get_center();
    std::array<double, CH_SPACEDIM> center2 = AH2->get_center();

    // set as center of the merged AH the center of the 2 AHs provided
    std::array<double, CH_SPACEDIM> center_merged = {
        (center1[0] + center2[0]) / 2., (center1[1] + center2[1]) / 2.,
        (center1[2] + center2[2]) / 2.};

    add_ah(a_bhamr, center_merged, initial_guess_merged, a_params, a_stats,
           a_coords);
    a_bhamr.m_apparent_horizons[num_ah]->set_merger_pair(ah1, ah2);
}

void AHFinder::solve(BHAMR &a_bhamr, double a_dt, double a_time,
                     double a_restart_time)
{
    CH_TIME("AHFinder::solve");

    // solve if it has never converged before OR if has already converged in the
    // past this means it will stop solving as soon as it stops converging but
    // keeps trying (with same initial guess) if never found one
    // (e.g. in the BinaryBH case, AH1 and AH2 stop solving once they
    // disappear, leaving the merged AH alone from then on)
    // (e.g. in the ScalarField case, only after a few time steps does an AH
    // form, and from then on it keeps solving forever)

    int num_ah = a_bhamr.m_apparent_horizons.size();
    auto &ahs = a_bhamr.m_apparent_horizons;

    // first solve for non-mergers
    for (int i = 0; i < num_ah; ++i)
    {
        if (ahs[i]->is_merger())
            continue; // solve non-mergers first
        if (!ahs[i]->stop_solving())
            ahs[i]->solve(a_dt, a_time, a_restart_time);
    }

    // now update the mergers and solve
    for (int i = 0; i < num_ah; ++i)
    {
        if (!ahs[i]->is_merger())
            continue; // solve mergers only

        // update center of merged if not yet converged, otherwise it does it by
        // itself in solve
        if (!ahs[i]->get_converged())
        {
            auto center1 = ahs[ahs[i]->get_merger_pair().first]->get_center();
            auto center2 = ahs[ahs[i]->get_merger_pair().second]->get_center();
            center1[0] = (center1[0] + center2[0]) / 2.;
            center1[1] = (center1[1] + center2[1]) / 2.;
            center1[2] = (center1[2] + center2[2]) / 2.;
            ahs[i]->set_center(center1);
        }

        // if (!ahs[i]->has_been_found() || ahs[i]->get_converged())
        ahs[i]->solve(a_dt, a_time, a_restart_time);
    }
}

/////////////////////////////////////////////////////////
// PETSc control methods
/////////////////////////////////////////////////////////

void AHFinder::set_num_ranks(int a_num_ranks)
{
    if (m_initialized)
        return;

#ifdef CH_MPI
    if (m_mpi_group != MPI_GROUP_NULL)
        MPI_Group_free(&m_mpi_group);

    if (m_mpi_comm != MPI_COMM_NULL)
        MPI_Comm_free(&m_mpi_comm);

    if (a_num_ranks > 0) // otherwise use the whole Chombo communicator
    {
        // don't make more ranks than there are processes
        int size;
        MPI_Comm_size(Chombo_MPI::comm, &size);
        a_num_ranks = std::min(a_num_ranks, size);

        int rank;
        MPI_Comm_rank(Chombo_MPI::comm, &rank);
        if (rank == 0)
            std::cout << "Using PETSc with " << a_num_ranks << " ranks"
                      << std::endl;

        MPI_Group MPI_GROUP_WORLD;
        MPI_Comm_group(Chombo_MPI::comm, &MPI_GROUP_WORLD);

        // (TF): for things like 'SmallDataIO', rank 0 of Chombo_MPI::comm
        // is the one that prints (not rank 0 from PETSC_COMM_WORLD) so pay
        // attention before taking rank 0 of Chombo out of PETSc. As of
        // 26/05/2019, I tested with 'int range[3] = { 1, a_num_ranks , 1 }'
        // and all the AH stuff was compatible with doing so
        int range[3] = {0, a_num_ranks - 1, 1};
        MPI_Group_range_incl(MPI_GROUP_WORLD, 1,
                             reinterpret_cast<int(*)[3]>(&range), &m_mpi_group);
        MPI_Group_free(&MPI_GROUP_WORLD);

        MPI_Comm_create(Chombo_MPI::comm, m_mpi_group, &m_mpi_comm);
    }
#endif
}

//! initialize PETSc and its MPI sub-communicator
PetscErrorCode AHFinder::initialize(int a_num_ranks, int argc, char *argv[])
{
    set_num_ranks(a_num_ranks);

    PetscErrorCode err = 0;
    if (!m_initialized)
    {
#ifdef CH_MPI
        if (m_mpi_comm != MPI_COMM_NULL)
            PETSC_COMM_WORLD = m_mpi_comm;
        else // use Chombo communicator if no 'set_num_ranks' was called (or
             // if called with 'a_num_ranks'<=0)
            PETSC_COMM_WORLD = Chombo_MPI::comm;
#endif

        if (AHFinder::is_rank_active())
        {
            if (argc > 0)
                err = PetscInitialize(&argc, &argv, NULL, NULL);
            else
                err = PetscInitializeNoArguments();
        }
        if (!err)
            m_initialized = true;
    }
    return err;
}

//! finalize PETSc
PetscErrorCode AHFinder::finalize()
{
    PetscErrorCode err = 0;
    if (m_initialized)
    {
        if (AHFinder::is_rank_active())
            err = PetscFinalize();

        if (!err)
            m_initialized = false;
    }
    return err;
}

//! true if part of PETSc MPI sub-communicator
bool AHFinder::is_rank_active()
{
#ifdef CH_MPI
    if (m_mpi_group == MPI_GROUP_NULL)
        return true;
    else
    {
        int rank;
        MPI_Group_rank(m_mpi_group, &rank);
        return rank != MPI_UNDEFINED;
    }
#else
    return true;
#endif
}

/////////////////////////////////////////////////////////
// initialize some "static" variables of AHFinder class
/////////////////////////////////////////////////////////

bool AHFinder::m_initialized = false;

#ifdef CH_MPI
MPI_Group AHFinder::m_mpi_group = MPI_GROUP_NULL;
MPI_Comm AHFinder::m_mpi_comm = MPI_COMM_NULL;
#endif

#endif
