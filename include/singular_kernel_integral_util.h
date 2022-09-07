#ifndef singular_kernel_integral_util_h
#define singular_kernel_integral_util_h

#include <deal.II/base/point.h>

using namespace dealii;

class singular_kernel_integral_util
{
public:
  ///
  /// Distance from point c to segment ab
  /// Page 130 in "Real time collision detection, Christer Ericson, Morgan Kaufmann, 2005".
  ///
  static double
  distance_squared_point_to_segment(const Point<3> &a, const Point<3> &b, const Point<3> &c);

  static double
  angle_between_two_vectors(const Tensor<1, 3> &a, const Tensor<1, 3> &b);

  ///
  /// Intersection point between line segments
  /// https://blogs.sas.com/content/iml/2018/07/09/intersection-line-segments.html
  ///
  static Tensor<1, 3>
  intersection_point_between_segments(const Tensor<1, 3> &p1, const Tensor<1, 3> &p2, const Tensor<1, 3> &q1, const Tensor<1, 3> &q2);

  ///
  /// Equation of external contour in polar coordinates
  ///     rho = rho(theta).
  /// see Figure 4 in Guiggiani.
  /// This version is based on law of sines.
  ///
  static double
  equation_of_external_contour_polar_coords_ver1(double theta_0, double theta, const Tensor<1, 3> &eta, const Tensor<1, 3> &v0, const Tensor<1, 3> &v1);

  ///
  /// Equation of external contour in polar coordinates
  ///     rho = rho(theta).
  /// see Figure 4 in Guiggiani.
  /// This version is based on intersection between line segments.
  ///
  static double
  equation_of_external_contour_polar_coords_ver2(double theta_0, double theta, const Tensor<1, 3> &eta, const Tensor<1, 3> &v0, const Tensor<1, 3> &v1);

  static void
  parameter_space_angles(const Point<3> &v0, const Point<3> &v1, const Point<3> &eta, double &theta_0, double &theta_1);

  static void
  expansion_functions_Fm1_Fm2(const Tensor<1, 3> &Jk0, const Tensor<1, 3> &Jk1, const Tensor<1, 3> &A, const Tensor<1, 3> &B, double q_A, double q_C, double N0, double N1, Tensor<1, 3> &FF_1, Tensor<1, 3> &FF_2);
};

#endif // singular_kernel_integral_util_h
