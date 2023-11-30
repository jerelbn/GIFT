/*
  This file is part of GIFT.

  GIFT is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  GIFT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with GIFT.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "GIFT/camera/GICamera.h"

namespace GIFT {

class CataCamera : public GICamera {
  /**
   * C. Mei, and P. Rives, Single View Point Omnidirectional Camera Calibration
   * from Planar Grids, ICRA 2007
   */
  protected:
    double m_focalLength;
    double m_xi;
    double m_k1;
    double m_k2;
    double m_p1;
    double m_p2;
    double m_gamma1;
    double m_gamma2;
    double m_u0;
    double m_v0;
    double m_inv_K11;
    double m_inv_K13;
    double m_inv_K22;
    double m_inv_K23;

    virtual Eigen::Vector2T projectPointEigen(const Eigen::Vector3T& point) const override;
    virtual Eigen::Vector3T undistortPointEigen(const Eigen::Vector2T& point) const override;

  public:
    CataCamera(int w, int h, double f, double xi, double k1, double k2, double p1, double p2, double gamma1, double gamma2, double u0, double v0);

    /** 
     * @brief Apply distortion to input point (from the normalised plane)
     *  
     * @param p_u undistorted coordinates of point on the normalised plane
     * @return the distorted point: p_d = p_u + d_u
     */
    void distortion(const Eigen::Vector2d& p_u, Eigen::Vector2d& d_u) const;
    void distortion(const Eigen::Vector2d& p_u, Eigen::Vector2d& d_u, Eigen::Matrix2d& J) const;

    virtual Eigen::Matrix<ftype, 2, 3> projectionJacobian(const Eigen::Vector3T& point) const override;
};

} // namespace GIFT
