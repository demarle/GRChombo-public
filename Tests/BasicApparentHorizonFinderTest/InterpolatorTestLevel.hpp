/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTERPOLATORTESTLEVEL_HPP_
#define INTERPOLATORTESTLEVEL_HPP_

#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FourthOrderFillPatch.H"
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "GRLevelData.hpp"
#include "InterpSource.hpp"
#include "LevelFluxRegister.H" //We don't actually use flux conservation but Chombo assumes we do
#include "LevelRK4.H"

class InterpolatorTestLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<InterpolatorTestLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual void initialData()
    {

        m_state_new.setVal(0);
        m_state_new.setVal(1, c_h11);
        m_state_new.setVal(1, c_h22);
        m_state_new.setVal(1, c_h33);
        const DisjointBoxLayout &level_domain = m_state_new.disjointBoxLayout();

        DataIterator dit = level_domain.dataIterator();
        double L = 16;
        double center = L / 2.;
        int N_GRID = 64;
        double mass = m_p.mass;

        for (dit.begin(); dit.ok(); ++dit)
        {
            FArrayBox &fab = m_state_new[dit];
            const Box &b = level_domain[dit];

            const IntVect &smallEnd = b.smallEnd();
            const IntVect &bigEnd = b.bigEnd();

            const int xmin = smallEnd[0];
            const int ymin = smallEnd[1];
            const int zmin = smallEnd[2];

            const int xmax = bigEnd[0];
            const int ymax = bigEnd[1];
            const int zmax = bigEnd[2];

            for (int iz = zmin - 3; iz <= zmax + 3; ++iz)
                for (int iy = ymin - 3; iy <= ymax + 3; ++iy)
                    for (int ix = xmin - 3; ix <= xmax + 3; ++ix)
                    {
                        const double zz = (iz + 0.5) * m_dx - center;
                        const double yy = (iy + 0.5) * m_dx - center;
                        const double xx = (ix + 0.5) * m_dx - center;
                        const double rr = sqrt(xx * xx + yy * yy + zz * zz);
                        const IntVect iv(ix, iy, iz);
                        fab(iv, c_chi) = pow((1 + mass / (2 * rr)), -4.0);
                    }
        }
    }

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
    {
    }

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state){};
};

#endif /* INTERPOLATORTESTLEVEL_HPP_ */
