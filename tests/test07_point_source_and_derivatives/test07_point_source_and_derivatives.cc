#include <cmath>
#include <iostream>

#include "laplace_kernel.h"

int
main(int argc, char **argv)
{
  double pi      = std::acos(-1.0);
  double four_pi = 4.0 * pi;

  double x0 = 0.0;
  double y0 = 0.0;
  double z0 = 0.0;

  double x = 1.0;
  double y = 0.0;
  double z = 0.0;

  double dx = x - x0;
  double dy = y - y0;
  double dz = z - z0;

  double dx2 = dx * dx;
  double dy2 = dy * dy;
  double dz2 = dz * dz;

  // Low speed aero dynamics equation 3.46-3.48:
  double denom = four_pi * std::pow(dx2 + dy2 + dz2, 5.0 / 2.0);
  double u     = -(dy2 + dz2 - 2 * dx2) / denom;
  double v     = 3.0 * dx * dy / denom;
  double w     = 3.0 * dx * dz / denom;

  std::cout << u << " " << v << " " << w << std::endl;


  const int dim = 3;
  dealii::Tensor<1, dim> direction;
  direction[0] = 0.0;
  direction[1] = 1.0;
  direction[2] = 0.0;
  dealii::Tensor<1, dim> normal;
  normal[0] = 0.0;
  normal[1] = 1.0;
  normal[2] = 0.0;
  double                 G0;
  dealii::Tensor<1, dim> G1;
  dealii::Tensor<1, dim> G2;
  LaplaceKernel::kernels(direction, normal, G0, G1, G2);

  std::cout << G2[0] << " " << G2[1] << " " << G2[2] << std::endl;
}
