#include "../include/singular_kernel_integral_util.h"

#include <deal.II/lac/full_matrix.h>

double
singular_kernel_integral_util::distance_squared_point_to_segment(const Point<3> &a, const Point<3> &b, const Point<3> &c)
{
  auto ab = b - a;
  auto ac = c - a;
  auto bc = c - b;

  auto e = ac * ab;

  if (e <= 0.0)
    return ac * ac;

  auto f = ab * ab;

  if (e >= f)
    return bc * bc;

  return ac * ac - e * e / f;
}

double
singular_kernel_integral_util::angle_between_two_vectors(const Tensor<1, 3> &a, const Tensor<1, 3> &b)
{
  return acos(a * b / (a.norm() * b.norm()));
}

///
/// Intersection point between line segments
/// https://blogs.sas.com/content/iml/2018/07/09/intersection-line-segments.html
///
Tensor<1, 3>
singular_kernel_integral_util::intersection_point_between_segments(const Tensor<1, 3> &p1, const Tensor<1, 3> &p2, const Tensor<1, 3> &q1, const Tensor<1, 3> &q2)
{
  auto p1p2 = p2 - p1;
  auto q2q1 = q1 - q2;
  auto p1q1 = q1 - p1;

  FullMatrix<double> A(3, 3);
  FullMatrix<double> Ainv(3, 3);
  A(0, 0) = p1p2[0];
  A(1, 0) = p1p2[1];
  A(0, 1) = q2q1[0];
  A(1, 1) = q2q1[1];
  A(2, 2) = 1.0;

  Ainv.invert(A);
  Tensor<2, 3> AAinv;
  Ainv.copy_to(AAinv);

  auto sol = AAinv * p1q1;

  auto s = sol[0];
  auto t = sol[1];

  auto r1 = (1.0 - s) * p1 + s * p2;
  auto r2 = (1.0 - t) * q1 + t * q2;

  return r1;
}

///
/// Equation of external contour in polar coordinates
///     rho = rho(theta).
/// see Figure 4 in Guiggiani.
/// This version is based on law of sines.
///
double
singular_kernel_integral_util::equation_of_external_contour_polar_coords_ver1(double theta_0, double theta, const Tensor<1, 3> &eta, const Tensor<1, 3> &v0, const Tensor<1, 3> &v1)
{
  auto r0  = v0 - eta;
  auto A   = theta - theta_0;
  auto B   = angle_between_two_vectors(v0 - v1, r0);
  auto C   = numbers::PI - A - B;
  auto c   = r0.norm();
  auto k   = c / sin(C);
  auto rho = k * sin(B);
  return rho;
}

///
/// Equation of external contour in polar coordinates
///     rho = rho(theta).
/// see Figure 4 in Guiggiani.
/// This version is based on intersection between line segments.
///
double
singular_kernel_integral_util::equation_of_external_contour_polar_coords_ver2(double theta_0, double theta, const Tensor<1, 3> &eta, const Tensor<1, 3> &v0, const Tensor<1, 3> &v1)
{
  auto r0  = v0 - eta;
  auto A   = theta - theta_0;
  auto tmp = eta + Point<3>(r0[0] * cos(A) - r0[1] * sin(A), r0[0] * sin(A) + r0[1] * cos(A), 0.0);
  auto p   = intersection_point_between_segments(eta, tmp, v0, v1);
  auto rho = (eta - p).norm();
  return rho;
}

void
singular_kernel_integral_util::parameter_space_angles(const Point<3> &v0, const Point<3> &v1, const Point<3> &eta, double &theta_0, double &theta_1)
{
  auto r0 = v0 - Point<3>(eta[0], eta[1], 0.0);
  auto r1 = v1 - Point<3>(eta[0], eta[1], 0.0);

  theta_0 = atan2(r0[1], r0[0]);
  if (theta_0 < 0)
    theta_0 += 2 * dealii::numbers::PI;

  theta_1 = atan2(r1[1], r1[0]);
  if (theta_1 < 0)
    theta_1 += 2 * dealii::numbers::PI;

  if (theta_1 < theta_0)
    theta_1 += 2 * dealii::numbers::PI;
}

void
singular_kernel_integral_util::expansion_functions_Fm1_Fm2(const Tensor<1, 3> &Jk0, const Tensor<1, 3> &Jk1, const Tensor<1, 3> &A, const Tensor<1, 3> &B, double q_A, double q_C, double N0, double N1, Tensor<1, 3> &FF_1, Tensor<1, 3> &FF_2)
{
  // (r,k  Jk) expansion
  // verified against (R/r)*inner_face_fe_values.JxW(p)*inner_uv_q_normals[p]
  double p_1 = (Jk0 * B + Jk1 * A) / q_A;

  // r,i expansion
  // verified against R/r
  Tensor<1, 3> d_0 = A / q_A;
  //              Tensor<1, 3> d_1 = B / q_A - A / pow(q_A, 3) * q_C;

  // r,i (r,k  Jk) expansion
  // verified against (R/r)*((R/r)*inner_face_fe_values.JxW(p)*inner_uv_q_normals[p])
  Tensor<1, 3> g_1 = d_0 * p_1;

  // 3*r,i (r,k  Jk) - Ji  expansion
  // verified against 3.0*(R/r)*((R/r)*inner_face_fe_values.JxW(p)*inner_uv_q_normals[p])-inner_uv_q_normals[p]*inner_face_fe_values.JxW(p)
  Tensor<1, 3> b_0 = -Jk0;
  Tensor<1, 3> b_1 = 3.0 * g_1 - Jk1;

  // 1/r^3  expansion
  // verified against 1/pow(r,3.0)
  double s_3 = 1.0 / pow(q_A, 3.0);
  double s_2 = -3.0 * q_C / pow(q_A, 5.0);

  const auto a_0 = b_0 * N0;
  const auto a_1 = b_1 * N0 + b_0 * N1;

  // -1/r^3 * (3*r,i (r,k  Jk) - Ji ) expansion
  FF_2 = -(s_3 * a_0);
  FF_1 = -(s_2 * a_0 + s_3 * a_1);
}