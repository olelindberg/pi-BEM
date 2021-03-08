
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
  single_layer(const Point<dim> &R)
  {
    switch (dim)
      {
        case 2:
          return (-std::log(R.norm()) / (2 * numbers::PI));

        case 3:
          return (1. / (R.norm() * 4 * numbers::PI));

        default:
          Assert(false, ExcInternalError());
          return 0.;
      }
  }



  template <int dim>
  Point<dim>
  double_layer(const Point<dim> &R)
  {
    switch (dim)
      {
        case 2:
          return R / (-2 * numbers::PI * R.square());
        case 3:
          return R / (-4 * numbers::PI * R.square() * R.norm());

        default:
          Assert(false, ExcInternalError());
          return Point<dim>();
      }
  }

  template <int dim>
  void
  kernels(const Tensor<1, dim> &R, Tensor<1, dim> &D, double &d)
  {
    double r  = R.norm();
    double r2 = r * r;
    switch (dim)
      {
        case 2:
          d = -std::log(r) / (2 * numbers::PI);
          D = R / (-2 * numbers::PI * r2);
          break;
        case 3:
          d = (1. / (r * 4 * numbers::PI));
          D = R / (-4 * numbers::PI * r2 * r);
          break;
        default:
          Assert(false, ExcInternalError());
      }
  }

  template <int dim>
  void
  double_body_kernel(const Tensor<1, dim> &x,const Tensor<1, dim> &x0, Tensor<1, dim> &D, double &d)
  {
//    switch (dim)
//      {
//        case 2:
//          break;
//        case 3:

          // Non reflected kernel:
          Tensor<1, dim> R = x0-x;
          kernels(R,D,d);

          // Reflected kernel:
          R     = x0;
          R[2] *= -1.0;
          R    -= x;

          double d1;
          Tensor<1, dim> D1;
          kernels(R,D1,d1);
          d += d1;
          D += D1;

//          break;
//        default:
//          Assert(false, ExcInternalError());
//      }

  }




} // namespace LaplaceKernel
