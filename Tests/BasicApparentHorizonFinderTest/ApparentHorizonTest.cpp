/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __    _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to LICENSE, in Chombo's root directory.
 */
#endif

// General includes:
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>

#include "parstream.H" //Gives us pout()
using std::endl;
#include "BHAMR.hpp"

#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "AMRInterpolator.hpp"
#include "ApparentHorizonTestLevel.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "UserVariables.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

int runApparentHorizonTest(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    GRParmParse pp(0, argv + argc, NULL, in_string.c_str());
    SimulationParameters sim_params(pp);

    BHAMR gr_amr;
    DefaultLevelFactory<ApparentHorizonTestLevel> ah_test_level_fact(
        gr_amr, sim_params);
    setupAMRObject(gr_amr, ah_test_level_fact);

    int status = 0;

#ifdef USE_AHFINDER
    AHFinder::params AH_params = {1, 40, 20, 1, 1, true, true, true};

    // Set up interpolator and PETSc subcommunicator when AH extraction is
    // active
    AMRInterpolator<Lagrange<4>> interpolator(
        gr_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    gr_amr.set_interpolator(&interpolator);

    AHFinder::add_ah(gr_amr, sim_params.center, sim_params.initial_guess,
                     AH_params);

    // get area from file to determine status
    auto stats = SmallDataIO::read("stats_AH1.dat");
    CH_assert(stats.size() > 0);
    double calculated_mass = stats[3][0];
    double error_perc =
        fabs(1. - calculated_mass / sim_params.kerr_params.mass) * 100;
    pout() << "error = " << error_perc << "%" << std::endl;
    status |= (error_perc > 0.1) ||
              std::isnan(error_perc); // accept 0.1% error in area calculation
#endif

    return status;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runApparentHorizonTest(argc, argv);

    if (status == 0)
    {
        std::cout << "BasicApparentHorizon test passed." << endl;
        pout() << "BasicApparentHorizon test passed." << endl;
    }
    else
    {
        std::cout << "BasicApparentHorizon test FAILED with return code "
                  << status << endl;
        pout() << "BasicApparentHorizon test FAILED with return code " << status
               << endl;
    }

    mainFinalize();
    return status;
}
