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

#include <opencv2/core/core.hpp>

#include "GIFT/camera/CataCamera.h"

using namespace GIFT;

CataCamera::CataCamera(int w, int h, double f, double xi, double k1, double k2, double p1, double p2, double gamma1, double gamma2, double u0, double v0)
 : m_focalLength(f)
 , m_xi(xi)
 , m_k1(k1)
 , m_k2(k2)
 , m_p1(p1)
 , m_p2(p2)
 , m_gamma1(gamma1)
 , m_gamma2(gamma2)
 , m_u0(u0)
 , m_v0(v0) {
  imageSize = cv::Size(w, h);
  m_inv_K11 = 1.0 / gamma1;
  m_inv_K13 = -u0 / gamma1;
  m_inv_K22 = 1.0 / gamma2;
  m_inv_K23 = -v0 / gamma2;
}

Eigen::Vector2T CataCamera::projectPointEigen(const Eigen::Vector3T& point) const {
  Eigen::Vector2d p_u, p_d;

  // Project points to the normalised plane
  double z = point(2) + m_xi * point.norm();
  p_u << point(0) / z, point(1) / z;

  // Apply distortion
  Eigen::Vector2d d_u;
  distortion(p_u, d_u);
  p_d = p_u + d_u;

  // Apply generalised projection matrix
  return Eigen::Vector2T(m_gamma1 * p_d(0) + m_u0,
                         m_gamma2 * p_d(1) + m_v0);
}

Eigen::Vector3T CataCamera::undistortPointEigen(const Eigen::Vector2T& point) const {
  // Lift points to normalised plane
  double mx_d = m_inv_K11 * point(0) + m_inv_K13;
  double my_d = m_inv_K22 * point(1) + m_inv_K23;

  // Recursive distortion model
  int n = 8;
  Eigen::Vector2d d_u;
  distortion(Eigen::Vector2d(mx_d, my_d), d_u);
  double mx_u = mx_d - d_u(0);
  double my_u = my_d - d_u(1);
  for (int i = 1; i < n; ++i) {
    distortion(Eigen::Vector2d(mx_u, my_u), d_u);
    mx_u = mx_d - d_u(0);
    my_u = my_d - d_u(1);
  }

  // Obtain a projective ray
  Eigen::Vector3T result;
  double rho2_d;
  if (m_xi == 1.0) {
    result << mx_u, my_u, (1.0 - mx_u * mx_u - my_u * my_u) / 2.0;
  } else {
    // Reuse variable
    rho2_d = mx_u * mx_u + my_u * my_u;
    result << mx_u, my_u, 1.0 - m_xi * (rho2_d + 1.0) / (m_xi + sqrt(1.0 + (1.0 - m_xi * m_xi) * rho2_d));
  }

  return result.normalized();
}

void CataCamera::distortion(const Eigen::Vector2d& p_u, Eigen::Vector2d& d_u) const {
  double mx2_u, my2_u, mxy_u, rho2_u, rad_dist_u;

  mx2_u = p_u(0) * p_u(0);
  my2_u = p_u(1) * p_u(1);
  mxy_u = p_u(0) * p_u(1);
  rho2_u = mx2_u + my2_u;
  rad_dist_u = m_k1 * rho2_u + m_k2 * rho2_u * rho2_u;

  d_u << p_u(0) * rad_dist_u + 2.0 * m_p1 * mxy_u + m_p2 * (rho2_u + 2.0 * mx2_u),
         p_u(1) * rad_dist_u + 2.0 * m_p2 * mxy_u + m_p1 * (rho2_u + 2.0 * my2_u);
}

void CataCamera::distortion(const Eigen::Vector2d& p_u, Eigen::Vector2d& d_u, Eigen::Matrix2d& J) const {
  double mx2_u, my2_u, mxy_u, rho2_u, rad_dist_u;

  mx2_u = p_u(0) * p_u(0);
  my2_u = p_u(1) * p_u(1);
  mxy_u = p_u(0) * p_u(1);
  rho2_u = mx2_u + my2_u;
  rad_dist_u = m_k1 * rho2_u + m_k2 * rho2_u * rho2_u;

  d_u << p_u(0) * rad_dist_u + 2.0 * m_p1 * mxy_u + m_p2 * (rho2_u + 2.0 * mx2_u),
         p_u(1) * rad_dist_u + 2.0 * m_p2 * mxy_u + m_p1 * (rho2_u + 2.0 * my2_u);

  double dxdmx = 1.0 + rad_dist_u + m_k1 * 2.0 * mx2_u + m_k2 * rho2_u * 4.0 * mx2_u + 2.0 * m_p1 * p_u(1) + 6.0 * m_p2 * p_u(0);
  double dydmx = m_k1 * 2.0 * p_u(0) * p_u(1) + m_k2 * 4.0 * rho2_u * p_u(0) * p_u(1) + m_p1 * 2.0 * p_u(0) + 2.0 * m_p2 * p_u(1);
  double dxdmy = dydmx;
  double dydmy = 1.0 + rad_dist_u + m_k1 * 2.0 * my2_u + m_k2 * rho2_u * 4.0 * my2_u + 6.0 * m_p1 * p_u(1) + 2.0 * m_p2 * p_u(0);

  J << dxdmx, dxdmy,
       dydmx, dydmy;
}

// TODO: fix and test this
Eigen::Matrix<ftype, 2, 3> CataCamera::projectionJacobian(const Eigen::Vector3T& point) const {
  // Numerical jacobian
  const ftype eps = 1e-5;
  Eigen::Matrix<ftype, 2, 3> J;
  for (int j = 0; j < 2; ++j) {
    const Eigen::Vector2T pp = projectPointEigen(point + eps * Eigen::Matrix3T::Identity().col(j));
    const Eigen::Vector2T pm = projectPointEigen(point - eps * Eigen::Matrix3T::Identity().col(j));
    J.col(j) = (pp - pm) / (2.0 * eps);
  }
  return J;
}
