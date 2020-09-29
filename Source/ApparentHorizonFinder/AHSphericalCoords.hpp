/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHSPHERICALCOORDS_HPP_
#define _AHSPHERICALCOORDS_HPP_

//! Class for coordinate system to use
//! Spherical coordinates
class AHSphericalCoords
{
  private:
    // center of sphere
    double m_center_x;
    double m_center_y;
    double m_center_z;

  public:
    AHSphericalCoords(const std::array<double, CH_SPACEDIM> &a_center)
    {
        set_center(a_center);
    }
    AHSphericalCoords(double a_center_x = 0, double a_center_y = 0,
                      double a_center_z = 0)
    {
        set_center(a_center_x, a_center_y, a_center_z);
    }

    void set_center(const std::array<double, CH_SPACEDIM> &a_center)
    {
        set_center(a_center[0], a_center[1], a_center[2]);
    }
    void set_center(double a_center_x, double a_center_y, double a_center_z)
    {
        m_center_x = a_center_x;
        m_center_y = a_center_y;
        m_center_z = a_center_z;
    }

    array<bool, CH_SPACEDIM - 1> is_periodic() const
    {
#if CH_SPACEDIM == 3
        const array<bool, 2> out = {false, true}; //{theta, phi}
#elif CH_SPACEDIM == 2
        const array<bool, 1> out = {false};       //{theta}
#endif
        return out;
    }

    //! upper bound of (theta,phi) coordinates
    array<double, CH_SPACEDIM - 1> get_domain_high() const
    {
#if CH_SPACEDIM == 3
        const array<double, 2> out = {M_PI, 2. * M_PI}; //{theta, phi}
#elif CH_SPACEDIM == 2
        const array<double, 1> out = {2. * M_PI}; //{theta}
#endif
        return out;
    }

    //! lower bound of (theta,phi) coordinates
    array<double, CH_SPACEDIM - 1> get_domain_low() const
    {
#if CH_SPACEDIM == 3
        const array<double, 2> out = {0., 0.}; //{theta, phi}
#elif CH_SPACEDIM == 2
        const array<double, 1> out = {0.};        //{theta}
#endif
        return out;
    }

    //! return labels of variables of this coordinate system (used for output
    //! files)
    std::vector<std::string> get_labels() const
    {
#if CH_SPACEDIM == 3
        return {"theta", "phi", "radius"}; //{theta, phi}
#elif CH_SPACEDIM == 2
        return {"theta", "radius"};               //{theta}
#endif
    }

    //! transform from spherical to cartesian
    //! u=theta, v=phi, f=sqrt_r
    array<double, 3> transform(double theta, double phi, double sqrt_r) const
    {
        CH_TIME("AHSphericalCoords::transform");

        const double r = sqrt_r * sqrt_r;

        const array<double, 3> out = {
            /* x = */ (r * sin(theta) * cos(phi)) + m_center_x,
            /* y = */ (r * sin(theta) * sin(phi)) + m_center_y,
            /* z = */ (r * cos(theta)) + m_center_z};

        return out;
    }

    //! transform derivatives from spherical to cartesian derivatives
    //! u=theta, v=phi, f=sqrt_r
    void transform_derivs(double theta, double phi, double sqrt_r,
                          AHGeometryData &out) const
    {
        CH_TIME("AHSphericalCoords::transform_derivs");

        const double r = sqrt_r * sqrt_r;
        const double d_sqrt_r = 0.5 / sqrt(r);
        const double dd_sqrt_r = -0.25 / (sqrt(r) * r);

        out.du[0] = (cos(theta) * cos(phi)) / r;
        out.du[1] = (cos(theta) * sin(phi)) / r;
        out.du[2] = -sin(theta) / r;

        out.dv[0] = (-sin(phi)) / (r * sin(theta));
        out.dv[1] = (cos(phi)) / (r * sin(theta));
        out.dv[2] = 0;

        double dfdx = sin(theta) * cos(phi);
        double dfdy = sin(theta) * sin(phi);
        double dfdz = cos(theta);

        out.ddu[0][0] = (cos(2 * theta) * cos(phi) * cos(phi) - cos(2 * phi)) /
                        (r * r * tan(theta));
        out.ddu[0][1] =
            ((cos(2 * theta) - 2) * sin(2 * phi)) / (2 * r * r * tan(theta));
        out.ddu[0][2] = -(cos(2 * theta) * cos(phi)) / (r * r);
        out.ddu[1][1] = (cos(2 * theta) * sin(phi) * sin(phi) + cos(2 * phi)) /
                        (r * r * tan(theta));
        out.ddu[1][2] = -(cos(2 * theta) * sin(phi)) / (r * r);
        out.ddu[2][2] = (sin(2 * theta)) / (r * r);
        out.ddu[1][0] = out.ddu[0][1];
        out.ddu[2][0] = out.ddu[0][2];
        out.ddu[2][1] = out.ddu[1][2];

        out.ddv[0][0] = sin(2 * phi) / (r * r * sin(theta) * sin(theta));
        out.ddv[0][1] = -cos(2 * phi) / (r * r * sin(theta) * sin(theta));
        out.ddv[0][2] = 0.;
        out.ddv[1][1] = -sin(2 * phi) / (r * r * sin(theta) * sin(theta));
        out.ddv[1][2] = 0.;
        out.ddv[2][2] = 0.;
        out.ddv[1][0] = out.ddv[0][1];
        out.ddv[2][0] = out.ddv[0][2];
        out.ddv[2][1] = out.ddv[1][2];

        double ddfdxdx = (cos(theta) * cos(theta) +
                          sin(theta) * sin(theta) * sin(phi) * sin(phi)) /
                         r;
        double ddfdxdy = -(sin(theta) * sin(theta) * sin(phi) * cos(phi)) / r;
        double ddfdxdz = -(sin(theta) * cos(theta) * cos(phi)) / r;
        double ddfdydy = (cos(theta) * cos(theta) +
                          sin(theta) * sin(theta) * cos(phi) * cos(phi)) /
                         r;
        double ddfdydz = -(sin(theta) * cos(theta) * sin(phi)) / r;
        double ddfdzdz = (sin(theta) * sin(theta)) / r;

        out.df[0] = d_sqrt_r * dfdx;
        out.df[1] = d_sqrt_r * dfdy;
        out.df[2] = d_sqrt_r * dfdz;

        out.ddf[0][0] = d_sqrt_r * ddfdxdx + dd_sqrt_r * dfdx * dfdx;
        out.ddf[0][1] = d_sqrt_r * ddfdxdy + dd_sqrt_r * dfdx * dfdy;
        out.ddf[0][2] = d_sqrt_r * ddfdxdz + dd_sqrt_r * dfdx * dfdz;
        out.ddf[1][1] = d_sqrt_r * ddfdydy + dd_sqrt_r * dfdy * dfdy;
        out.ddf[1][2] = d_sqrt_r * ddfdydz + dd_sqrt_r * dfdy * dfdz;
        out.ddf[2][2] = d_sqrt_r * ddfdzdz + dd_sqrt_r * dfdz * dfdz;
        out.ddf[1][0] = out.ddf[0][1];
        out.ddf[2][0] = out.ddf[0][2];
        out.ddf[2][1] = out.ddf[1][2];

        out.dxdu[0] = r * cos(theta) * cos(phi);
        out.dxdu[1] = r * cos(theta) * sin(phi);
        out.dxdu[2] = -r * sin(theta);

        out.dxdv[0] = -r * sin(theta) * sin(phi);
        out.dxdv[1] = r * sin(theta) * cos(phi);
        out.dxdv[2] = 0;

        out.dxdf[0] = sin(theta) * cos(phi) / d_sqrt_r;
        out.dxdf[1] = sin(theta) * sin(phi) / d_sqrt_r;
        out.dxdf[2] = cos(theta) / d_sqrt_r;
    }
};

#endif /* _AHSPHERICALCOORDS_HPP_ */
