/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHFINDER_HPP_
#define _AHFINDER_HPP_

// Chombo includes
#include "SPMD.H" //Chombo_MPI::comm

// Other includes:

// included so that user can re-define default AHSurfaceGeometry or AHFunction
#include "UserVariables.hpp"

// define SurfaceGeometry of AHFinder
#ifndef AHSurfaceGeometry
#include "AHSphericalGeometry.hpp"
#define AHSurfaceGeometry AHSphericalGeometry
#endif

// default to expansion
#ifndef AHFunction
#include "AHFunctions.hpp"
#define AHFunction ExpansionFunction
#endif

#ifdef CH_MPI
#include <mpi.h>
#endif

#include <petsc.h>

#include "AHData.hpp"
#include "AHDeriv.hpp"
#include "AHGeometryData.hpp"
#include "AMRInterpolator.hpp"
#include "BoundaryConditions.hpp"
#include "ChomboParameters.hpp"
#include "Lagrange.hpp"
#include "VariableType.hpp"

// Chombo namespace
#include "UsingNamespace.H"

// declare classes, define them later
template <class SurfaceGeometry, class AHFunction> class AHInterpolation;
template <class SurfaceGeometry, class AHFunction> class ApparentHorizon;

/*
TF:
A general note on the radius and the value of chi on the AH surface.
For a Kerr BH initial data with dimensionless spin 's' and mass 'M':
  r_AH   = M / 4 (1 + sqrt{1 - s^2})
  chi_AH = [(1-s^2)^(1/6) / 16 * ((1 + sqrt{1-s^2})/2)^{2/3} for theta=0,
            (1-s^2)^{1/6} / 16 * ((1 + sqrt{1-s^2})/2)^{1/3} for theta=pi/2]

After puncture gauge settled:
  r_AH   ~ M * (0.275 + 0.76 * sqrt{1-s^2})
  chi_AH ~ 0.26 * sqrt{1-s^2}

Take it with a grain of salt, but these results are from a numerical fit I did
for numerical Kerr BH runs, approximately settled at around t~50M, from spin=0
to spin=0.99. You can use it generically for other things, like realizing that
for spin=0 you can plot in VisIt the contour chi=0.26 to get the AH, or that for
spin=0 the AH goes from r=M/2 to r~M

Binaries:
After some inspiral cases I've been through, I realized that merger happened
when the distance between the punctures was about (roughtly!) ~[0.5,2]*(M1+M2)
and with a radius of ~[0.7,1.1]*(M1+M2) (wobbly of course, and bigger values for
smaller final spin). It makes sense to start looking for the merger before this,
but remember that for PETSc diverging is significantly slower than converging,
so the closer the better. Hence, as default, start looking for a merger when the
separation is about ~<2*(M1+M2) and with a initial guess of ~2*(M1+M2) (twice as
big as what we'll probably get to make sure it converges to the outer AH; even
then it might converge to the inner one, but it is much easer for the AH to
converge if the initial guess is bigger than the actual radius, but the opposite
is more senstitive).
*/

//! Class to manage AHs and its mergers + control PETSc MPI sub-communicator
class AHFinder
{
  public:
    struct params
    {
        int num_ranks; //!< number of ranks for PETSc sub-communicator (default
                       //!< 0, which is 'all')
        int num_points_u; //!< number of points for 2D coordinate grid
#if CH_SPACEDIM == 3
        int num_points_v; //!< number of points for 2D coordinate grid
#endif
        int solve_interval; //!< same as checkpoint_interval, for
                            //!< ApparentHorizon::solve (default 1)
        int print_interval; //!< same as solve_interval, but for prints (default
                            //!< 1)
        bool track_center;  //!< whether or not to update the center
                            //!< (set to false if you know it won't move)
                            //!< (default true)
        bool predict_origin; //!< whether or not to estimate where the next
                             //!< origin will be at (default = track_center)

        int level_to_run; // if negative, it will count backwards (e.g. -1 is
                          // 'last level') (default 0)

        double start_time;   //!< time after which AH can start (default 0.)
                             //!< Useful for ScalarField collapse
        double give_up_time; //!< stop if at this time nothing was found
                             //!< (<0 to never, which is default)
                             //!< Useful for ScalarField collapse

        //! mergers will be searched when distance between 'parent' BHs is
        //! distance < merger_search_factor * 4. * (AH_initial_guess_1 +
        //! AH_initial_guess_2) should be roughly '2M=2(m1+m2)' for initial
        //! guess at m/2 (set to non-positive to 'always search')
        double merger_search_factor; // see note above (default is 1)
        //! initial guess for merger is 'merger_pre_factor * 4. *
        //! (AH_initial_guess_1 + AH_initial_guess_2)'
        //! set to somethig bigger to avoid finding the inner AH
        double merger_pre_factor; // see note above (default to 1.)

        bool allow_re_attempt; //!< re-attempt with initial guess if
                               //!< previous convergence failed (default false)
        int max_fails_after_lost; //!< number of time steps to try again after
                                  //!< (-1 to never) the AH was lost
                                  //!< (default is 0)

        int verbose; //!< print information during execution (default is 1)

        bool print_geometry_data; //!< print metric and extrinsic
                                  //!< curvature of the AHs (default false)

        bool re_solve_at_restart; //!< whether or not to re-run even if AH
                                  //!< already exists (useful in order to be
                                  //!< able to provide an initial guess and
                                  //!< re-run the AHFinder on that time step)
                                  //!< (default false)

        bool stop_if_max_fails; //! breaks the run if AH doesn't converge
                                //! 'max_fails_after_lost' times or if
                                //! 'give_up_time' is reached without
                                //! convergence (default is 'false')

        std::map<std::string, std::tuple<int, VariableType, int>>
            extra_vars;     //! extra vars to write to coordinates file (<enum,
                            //! evolution or diagnostic, int for local|d1|d2>)
        int num_extra_vars; // total number of extra vars (!=extra_vars.size()
                            // as derivative count for multiple vars)

        int extra_contain_diagnostic; // not a parameter (set internally);
                                      // counts how many

        void read_params(GRParmParse &pp, const ChomboParameters &a_p);
    };

    enum verbose_level
    {
        NONE,
        MIN,  // minimal
        SOME, // a bit technical
        MORE  // debug
    };

  private:
    //! if this AH is supposed to track the formation of a merger
    //! this pair indicates the indices of the 2 AHs in the vector
    //! 'm_apparent_horizons'
    std::vector<std::pair<int, int>> m_merger_pairs;
    AMRInterpolator<Lagrange<4>> *m_interpolator; //!< The interpolator pointer

    std::vector<ApparentHorizon<AHSurfaceGeometry, AHFunction> *>
        m_apparent_horizons; //!< public in case user wants to solve by himself

  public:
    AHFinder(){};
    ~AHFinder();

    ALWAYS_INLINE ApparentHorizon<AHSurfaceGeometry, AHFunction> *
    get(unsigned AH_i)
    {
        CH_assert(AH_i < m_apparent_horizons.size());
        return m_apparent_horizons[AH_i];
    }

    //! returns the index of the AH in m_apparent_horizons
    int
    add_ah(const AHSurfaceGeometry &a_coord_system,
           double a_initial_guess, //!< Initial guess for radius (or whatever
                                   //!< coordinate you're solving for)
           const params &a_params, //!< set of AH parameters
           const std::string &a_stats =
               "stats_AH", //!< name for stats file with
                           //!< area, spin and AH origin/center
           const std::string &a_coords =
               "coords_AH",             //!< name for coords file with AH
                                        //!< coordinates at each time step)
           bool solve_first_step = true //!< whether or not to solve if t=0
    );

    //! returns the index of the AH in m_apparent_horizons
    //! allows for personalized optimizer that finds zero of function
    //! 'AHFunction' (that can have some ::params)
    int add_ah(const AHSurfaceGeometry &a_coord_system, double a_initial_guess,
               const params &a_params,
               const typename AHFunction::params &a_func_params,
               const std::string &a_stats = "stats",
               const std::string &a_coords = "coords",
               bool solve_first_step = true);

    // returns the index of the AH in m_apparent_horizons
    int add_ah_merger(int ah1, int ah2, const params &a_params);

    //! Find AH; Calculate area and spin; Update center; Print outputs
    void solve(double a_dt, double a_time, double a_restart_time);

    bool need_diagnostics(double a_dt, double a_time)
        const; //!< is any AH printing diagnostics? Useful if Diagnostics need
               //!< to be calculated

    ALWAYS_INLINE void
    set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }

    void
    set_origins(const std::vector<std::array<double, CH_SPACEDIM>> &origins,
                bool includes_mergers = false);

  private:
    //! returns false if 'parent' AHs are too far
    //! sets the initial guess and the origin for the merger
    bool solve_merger(int ah1, int ah2, double &initial_guess_merger,
                      std::array<double, CH_SPACEDIM> &origin_merged);

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    ///////////////// PETSc control methods /////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////

    static bool m_initialized; //!< is initialized?

#ifdef CH_MPI
    static MPI_Group m_mpi_group; //!< set of MPI ranks for PETSc
    static MPI_Comm m_mpi_comm;   //!< MPI sub-communicator
#endif

    //! define number of ranks of PETSc sub-communicator
    static void set_num_ranks(int a_num_ranks);

    //! initialize PETSc and its MPI sub-communicator
    static PetscErrorCode PETSc_initialize(int a_num_ranks);

    //! finalize PETSc
    static PetscErrorCode PETSc_finalize();

  public:
    //! true if part of PETSc MPI sub-communicator
    static bool is_rank_active();
}; // namespace AHFinder

#endif /* _AHFINDER_HPP_ */
