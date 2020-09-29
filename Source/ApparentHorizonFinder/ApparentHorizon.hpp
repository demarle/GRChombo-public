/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _APPARENTHORIZON_HPP_
#define _APPARENTHORIZON_HPP_

#include <petsc.h>
#include <petscviewerhdf5.h>

#include "AHFinder.hpp"

#define DWIDTH 5
#define DDWIDTH 6

template <typename AHGeom> class ApparentHorizon
{
  public:
    ApparentHorizon(
        AHGeom &geom, //!< Geometry class to exchange data
        const std::array<double, CH_SPACEDIM>
            &a_center,          //!< initial guess for center of AH
        double a_initial_guess, //!< Initial guess for radius (or whatever
                                //!< coordinate you're solving for)
        const AHFinder::params &a_params, //!< set of AH parameters
        const std::string &a_stats =
            "stats.out", //!< name for output file with area, spin and AH center
        const std::string &a_coords =
            "coords_" //!< name for output file with AH coordinates at each time
                      //!< step
    );
    ~ApparentHorizon();

    void solve(double a_dt, double a_time, double a_restart_time,
               bool a_first_step = false);
    bool get_converged()
        const; //!< PETSc didn't converge last time solve() was called
    bool has_been_found() const; //!< PETSc never converged so far
    bool stop_solving() const;
    std::array<double, CH_SPACEDIM> get_center() const;
    void set_center(const std::array<double, CH_SPACEDIM> &a_center);

    double get_max_F() const;
    double get_min_F() const;

    double get_initial_guess();

    //! set/get/know whether or not this AH comes from a merger
    //! (useful for AHFinder::solve)
    void set_merger_pair(int, int);
    std::pair<int, int> get_merger_pair();
    bool is_merger();

    // methods
  private:
    void
    restart(int a_int_step,
            double a_current_time); //!< restart AH, updating the coordinates
                                    //!< based on the last output file and the
                                    //!< center based on the
    double calculate_area();        //!< calculate AH area, ONLY FOR 3D
    double
    calculate_coord_area(); //!< calculate AH coordinate area, ONLY FOR 3D
    double calculate_spin(
        double a_area); //!< calculate spin with 'z' direction, ONLY FOR 3D
    std::array<double, CH_SPACEDIM>
    calculate_center(); //!< update location of center by calculating the
                        //!< centroid of the AH
    //! Calculates the global maximum of F and stores it in m_max_F
    void calculate_minmax_F();

    void set_initial_guess();

    bool do_solve(double a_dt, double a_time)
        const; //!< decide (based times passed to 'solve') whether or not
               //!< to print (uses AHFinder::params::solve_interval)
    bool do_print(double a_dt, double a_time)
        const; //!< decide when to print (only AHFinder::params::print_interval
               //!< out of all 'solve's)

    void write_coords_to_file(double a_dt, double a_time, double a_restart_time)
        const;                //!< write coords in (m_u, m_v, m_F) to 'filename'
    void check_convergence(); //!< check if PETSc has converged and share that
                              //!< info across all Chombo processes

    void initialise_PETSc(); //!< initialise automatically done in constructor
    void finalise_PETSc();   //!< finalise automatically done in destructor

    // variables
  private:
    const AHFinder::params &m_params; //!< set of AH parameters

    const std::string m_stats, m_coords; //!< base names for output files

    bool m_converged; //!< flag saying if PETSc has converged or not (read using
                      //!< 'get_converged()')
    bool m_has_been_found; //!< flag saying if an horizon has ever been found
                           //!< (read using 'has_been_found()')
    int m_num_failed_convergences; //!< the number of failed consecutive
                                   //!< convergences
    std::array<double, CH_SPACEDIM> m_center;     //!< center of AH
    std::array<double, CH_SPACEDIM> m_center_old; //!< the old center of the AH
    double m_initial_guess; //!< initial guess for AH (saved so that it can be
                            //!< re-used when atempting to solve again)

    //! if this AH is supposed to track the formation of a merger
    //! this pair indicates the indices of the 2 AHs
    //! (according to the vector in BHAMR::m_apparent_horizons)
    std::pair<int, int> m_merger_pair;

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////// PETSc stuff below ///////////////////////////
    /////////////////////////////////////////////////////////////////////////

#if CH_SPACEDIM == 3
    typedef PetscScalar **dmda_arr_t;
#elif CH_SPACEDIM == 2
    typedef PetscScalar *dmda_arr_t;
#endif

    //! Geometries of the AH
    //! 'm_geom_plus' and 'm_geom_minus' are used to calculate the
    //! jacobian of the expansion using a 'delta' numerical differentiation
    AHGeom m_geom;
    AHGeom m_geom_plus;
    AHGeom m_geom_minus;

    //!< used to compute jacobian of expansion (numerical differentiation)
    static constexpr double eps = 1e-7;

    const bool m_periodic_u; //!< is 'u' periodic?
    int m_num_global_u;      //!< total number of grid points in 'u' coordinate
    double m_du;             //!< physical 'delta' in 'u' coordinate

    const bool m_periodic_v; //!< is 'v' periodic?
    int m_num_global_v;      //!< total number of grid points in 'u' coordinate
    double m_dv;             //!< physical 'delta' in 'v' coordinate

    //! vectors to store and manipulate 'u', 'v' and 'F'
    //! internally, in interaction with the 'AHGeom's
    std::vector<double> m_u;
    std::vector<double> m_v;
    std::vector<double> m_F;

    //! stores the global maximum and minimum of F - calculate with
    //! calculate_minmax_F
    double m_max_F;
    double m_min_F;

    DM m_dmda;

    //! minimums and maximums of coordinates 'u' and 'v'
    //! of the PETSc grid specific to the current rank
    int m_umin;
    int m_umax;

    int m_vmin;
    int m_vmax;

    //! number of points in 'u' and 'v' direction
    //! (m_nu = m_umax - m_umin)
    int m_nu;
    int m_nv;

    //! Scalable Nonlinear Equations Solvers
    SNES m_snes;

    Vec m_snes_soln;
    Vec m_snes_rhs;
    Mat m_snes_jac;

    Vec m_snes_guu;
    Vec m_snes_gvv;
    Vec m_snes_guv;
    Vec m_snes_gww;
    Vec m_snes_lapse;

    //! class to store 1st and 2nd derivatives of 'F'
    struct Deriv
    {
        double duF;
        double dvF;
        double duduF;
        double dvdvF;
        double dudvF;

        int du_stencil_start;
        int dv_stencil_start;
        int dudu_stencil_start;
        int dvdv_stencil_start;

        double du_weights[DWIDTH];
        double dv_weights[DWIDTH];
        double dudu_weights[DDWIDTH];
        double dvdv_weights[DDWIDTH];

        Deriv()
        {
            // force all (double) elements of Deriv to be 0
            memset(this, 0, sizeof(Deriv));
        }
    };

    //! function to calculate 1st and 2nd derivatives of 'in'
    //! (tipically corresponds to our 'f' function)
    //! in the 'u' and 'v' directions
    Deriv diff(dmda_arr_t in, int u, int v);

    //! private functions used to compute the RHS (the expansion) and it's
    //! jacobian
    void form_function(Vec F, Vec Rhs);
    void form_jacobian(Vec F, Mat J);

    //! functions used by PETSc based on 'form_function' and 'form_jacobian'
    static PetscErrorCode Petsc_form_function(SNES snes, Vec F, Vec Rhs,
                                              void *ptr)
    {
        ApparentHorizon<AHGeom> &ah =
            *reinterpret_cast<ApparentHorizon<AHGeom> *>(ptr);
        CH_assert(ah.m_snes == snes);
        ah.form_function(F, Rhs);
        return 0;
    }

    static PetscErrorCode
#if PETSC_VERSION_GE(3, 5, 0)
    Petsc_form_jacobian(SNES snes, Vec F, Mat Amat, Mat Pmat, void *ptr)
#else
    Petsc_form_jacobian(SNES snes, Vec F, Mat *Amat, Mat *Pmat,
                        MatStructure *flag, void *ptr)
#endif
    {
        ApparentHorizon<AHGeom> &ah =
            *reinterpret_cast<ApparentHorizon<AHGeom> *>(ptr);

#if PETSC_VERSION_GE(3, 5, 0)
        CH_assert(ah.m_snes == snes && Amat == Pmat);
        ah.form_jacobian(F, Amat);
#else
        CH_assert(ah.m_snes == snes && *Amat == *Pmat);
        ah.form_jacobian(F, *Amat);
#endif

        return 0;
    }

    static PetscErrorCode Petsc_SNES_monitor(SNES snes, PetscInt its,
                                             PetscReal norm, void *ptr)
    {
        ApparentHorizon<AHGeom> &ah =
            *reinterpret_cast<ApparentHorizon<AHGeom> *>(ptr);
        CH_assert(ah.m_snes == snes);
        return 0;
    }
};

#include "ApparentHorizon.impl.hpp"

#endif /* _APPARENTHORIZON_HPP_ */
