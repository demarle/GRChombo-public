/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "MatterConstraints.hpp"

// For tag cells
#include "ChiAndPhiExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "AMRReductions.hpp"
#include "ChiRelaxation.hpp"
#include "ComputePack.hpp"
#include "EMTensor.hpp"
#include "MatterEnergyFlux.hpp"
#include "MatterEnergyFluxExtraction.hpp"
#include "Potential.hpp"
#include "ScalarBubble.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"
#include "SetValueVolume.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero then initial conditions for scalar field -
    // here a bubble
    BoxLoops::loop(make_compute_pack(SetValue(0.0),
                                     ScalarBubble(m_p.initial_params, m_dx)),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    BoxLoops::loop(SetValue(0.0), m_state_diagnostics, m_state_diagnostics,
                   INCLUDE_GHOST_CELLS);
}

// Things to do before outputting a checkpoint file
void ScalarFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        make_compute_pack(MatterConstraints<ScalarFieldWithPotential>(
                              scalar_field, m_dx, m_p.G_Newton),
                          MatterEnergyFlux<ScalarFieldWithPotential>(
                              scalar_field, m_p.extraction_params.center, m_dx,
                              c_energy_flux)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{

    // Relaxation function for chi - this will eventually be done separately
    // with hdf5 as input
    if (m_time < m_p.relaxtime)
    {
        // Calculate chi relaxation right hand side
        // Note this assumes conformal chi and Mom constraint trivially
        // satisfied  No evolution in other variables, which are assumed to
        // satisfy constraints per initial conditions
        Potential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        ChiRelaxation<ScalarFieldWithPotential> relaxation(
            scalar_field, m_dx, m_p.relaxspeed, m_p.G_Newton);
        SetValue set_other_values_zero(0.0, Interval(c_h11, NUM_VARS - 1));
        auto compute_pack1 =
            make_compute_pack(relaxation, set_other_values_zero);
        BoxLoops::loop(compute_pack1, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else
    {

        // Enforce trace free A_ij and positive chi and alpha
        BoxLoops::loop(
            make_compute_pack(TraceARemoval(),
                              PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
            a_soln, a_soln, INCLUDE_GHOST_CELLS);

        // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
        Potential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        MatterCCZ4<ScalarFieldWithPotential> my_ccz4_matter(
            scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
            m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(ChiAndPhiExtractionTaggingCriterion(
                       m_dx, m_p.regrid_threshold_chi, m_p.regrid_threshold_phi,
                       m_p.extraction_params, m_level, m_p.activate_extraction),
                   current_state, tagging_criterion);
}

void ScalarFieldLevel::specificPostTimeStep()
{
    CH_TIME("ScalarFieldLevel::specificPostTimeStep");
    bool first_step = (m_dt == m_time);
    if (m_p.activate_extraction)
    {
        CH_TIME("energy_flux_extraction");
        fillAllGhosts();
        Potential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        BoxLoops::loop(
            make_compute_pack(
                EMTensor<ScalarFieldWithPotential>(scalar_field, m_dx, c_rho),
                MatterEnergyFlux<ScalarFieldWithPotential>(
                    scalar_field, m_p.extraction_params.center, m_dx,
                    c_energy_flux)),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            MatterEnergyFluxExtraction energy_flux_extraction(
                m_p.extraction_params, m_dt, m_time, first_step, m_restart_time,
                c_energy_flux);
            m_gr_amr.m_interpolator->refresh();
            energy_flux_extraction.execute_query(m_gr_amr.m_interpolator);
        }

        if (m_level == 0)
        {
            CH_TIME("calculate_energy");
            // Calculate the total energy within in each extraction sphere.
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            std::vector<GRAMRLevel *> all_level_ptrs =
                m_gr_amr.get_gramrlevels();
            std::vector<double> total_energy(
                m_p.extraction_params.num_extraction_radii);

            // start with largest radius first as we're going to zero the
            // grid variable outside each radius before computing the sum
            // note this assumes the extraction_radii are in ascending order
            // which is not enforced
            for (int iradius = m_p.extraction_params.num_extraction_radii - 1;
                 iradius >= 0; --iradius)
            {
                for (auto level_ptr : all_level_ptrs)
                {
                    ScalarFieldLevel *sf_level_ptr =
                        dynamic_cast<ScalarFieldLevel *>(level_ptr);
                    if (sf_level_ptr == nullptr)
                    {
                        break;
                    }

                    BoxLoops::loop(
                        SetValueVolume<SphericalVolumeComplement>(
                            0.0, Interval(c_rho, c_rho),
                            SphericalVolumeComplement(
                                m_p.extraction_params.extraction_radii[iradius],
                                m_p.extraction_params.center,
                                sf_level_ptr->m_dx)),
                        sf_level_ptr->m_state_diagnostics,
                        sf_level_ptr->m_state_diagnostics, EXCLUDE_GHOST_CELLS);
                }
                total_energy[iradius] = amr_reductions.sum(c_rho);
            }

            SmallDataIO energy_file("energy", m_dt, m_time, m_restart_time,
                                    SmallDataIO::APPEND, first_step);
            energy_file.remove_duplicate_time_data();
            if (first_step)
            {
                std::vector<std::string> header_strings(
                    m_p.extraction_params.num_extraction_radii);
                for (int iradius = 0;
                     iradius < m_p.extraction_params.num_extraction_radii;
                     ++iradius)
                {
                    header_strings[iradius] =
                        "r < " +
                        std::to_string(
                            m_p.extraction_params.extraction_radii[iradius]);
                }
                energy_file.write_header_line(header_strings);
            }
            energy_file.write_time_data_line(total_energy);
        }
    }
}
