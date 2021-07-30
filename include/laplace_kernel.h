
#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>

// And here are a few C++ standard header
// files that we will need:
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>


namespace LaplaceKernel
{
  template <int dim>
  double
  single_layer(const dealii::Point<dim> &R)
  {
    switch (dim)
      {
        case 2:
          return (-std::log(R.norm()) / (2 * dealii::numbers::PI));

        case 3:
          return (1. / (R.norm() * 4 * dealii::numbers::PI));

        default:
          Assert(false, ExcInternalError());
          return 0.;
      }
  }



  template <int dim>
  dealii::Point<dim>
  double_layer(const dealii::Point<dim> &R)
  {
    switch (dim)
      {
        case 2:
          return R / (-2 * dealii::numbers::PI * R.square());
        case 3:
          return R / (-4 * dealii::numbers::PI * R.square() * R.norm());

        default:
          Assert(false, ExcInternalError());
          return dealii::Point<dim>();
      }
  }

  template <int dim>
  void
  kernels(const dealii::Tensor<1, dim> &R, dealii::Tensor<1, dim> &D, double &d)
  {
    double r  = R.norm();
    double r2 = r * r;
    switch (dim)
      {
        case 2:
          d = -std::log(r) / (2 * dealii::numbers::PI);
          D = R / (-2 * dealii::numbers::PI * r2);
          break;
        case 3:
          d = (1. / (r * 4 * dealii::numbers::PI));
          D = R / (-4 * dealii::numbers::PI * r2 * r);
          break;
        default:
          Assert(false, ExcInternalError());
      }
  }

  template <int dim>
  void
  kernels(
    const dealii::Tensor<1,dim>& direction, 
    const dealii::Tensor<1,dim>& normal, 
    double& G0, 
    dealii::Tensor<1, dim>& G1, 
    dealii::Tensor<1, dim>& G2)
  {

    double two_pi  = 2.0*dealii::numbers::PI;
    double four_pi = 4.0*dealii::numbers::PI;

    double r  = direction.norm();
    double r2 = r * r;
    double r3 = r * r2;

    double dx = direction[0];
    double dy = direction[1];
    double dz = direction[2];

    double nx = normal[0];
    double ny = normal[1];
    double nz = normal[2];

    double dx2 = dx*dx;
    double dy2 = dy*dy;
    double dz2 = dz*dz;

    double rx  = dx/r;
    double ry  = dy/r;
    double rz  = dz/r;
    double rx2  = rx*rx;
    double ry2  = ry*ry;
    double rz2  = rz*rz;

    double rxx = (-dx2/r3 + 1.0/r);
    double ryy = (-dy2/r3 + 1.0/r);
    double rzz = (-dz2/r3 + 1.0/r);
    double rxy = -dx*dy/r3;
    double rxz = -dx*dz/r3;
    double ryz = -dy*dz/r3;
 
    G0    = 1.0/(four_pi*r);

    G1[0] = -rx/(four_pi*r2);
    G1[1] = -ry/(four_pi*r2);
    G1[2] = -rz/(four_pi*r2);
    
    G2[0] = -nx*rxx/(four_pi*r2) + nx*rx2/(two_pi*r3)   - ny*rxy/(four_pi*r2) + ny*rx*ry/(two_pi*r3) - nz*rxz/(four_pi*r2) + nz*rx*rz/(two_pi*r3);
    G2[1] = -nx*rxy/(four_pi*r2) + nx*rx*ry/(two_pi*r3) - ny*ryy/(four_pi*r2) + ny*ry2/(two_pi*r3)   - nz*ryz/(four_pi*r2) + nz*ry*rz/(two_pi*r3);
    G2[2] = -nx*rxz/(four_pi*r2) + nx*rx*rz/(two_pi*r3) - ny*ryz/(four_pi*r2) + ny*ry*rz/(two_pi*r3) - nz*rzz/(four_pi*r2) + nz*rz2/(two_pi*r3);

  }
} // namespace LaplaceKernel
