#ifndef singular_kernel_integral_h
#define singular_kernel_integral_h

// The most fundamental class in the library is the Triangulation class, which
// is declared here:
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria.h>
// And this is the file in which the functions are declared that create grids:
#include <deal.II/grid/grid_generator.h>

// This file contains the description of the Lagrange interpolation finite
// element:
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>


// And this file is needed for the creation of sparsity patterns of sparse
// matrices, as shown in previous examples:
#include <deal.II/dofs/dof_tools.h>

// The next two files are needed for assembling the matrix using quadrature on
// each cell. The classes declared in them will be explained below:
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_manifold.h>
// We need the following two includes for loops over cells and/or faces:
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
// Here are some functions to generate standard grids:
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>

// Output of grids in various graphics formats:
#include <deal.II/grid/grid_out.h>

// This is needed for C++ output:
#include <fstream>
#include <iostream>
// And this for the declarations of the `std::sqrt` and `std::fabs` functions:
#include <cmath>

// The final step in importing deal.II is this: All deal.II functions and
// classes are in a namespace <code>dealii</code>, to make sure they don't
// clash with symbols from other libraries you may want to use in conjunction
// with deal.II. One could use these functions and classes by prefixing every
// use of these names by <code>dealii::</code>, but that would quickly become
// cumbersome and annoying. Rather, we simply import the entire deal.II
// namespace for general use:
using namespace dealii;

///
/// The singular kernel integral is implemented based on [1] and [2].
///
/// References:
/// [1] M. Guiggiani, "Formulation and numerical treatment of boundary integral
/// equations with hypersingular kernels", 1998. [2] V. Mantic and F. Paris,
/// "Existence and evalution of the two free trems in the hypersingular boundary
/// integral equation of potential theory",
///     Engineering analysis with Boundary Elements 16 (1995) 253-260.
///
template <int dim>
class SingularKernelIntegral
{
public:
  SingularKernelIntegral(double in_rho_quadrature_order, double in_theta_quadrature_order, FiniteElement<dim - 1, dim> &in_fe, Mapping<dim - 1, dim> &in_mapping);

  Tensor<1, dim>
  evaluate_free_term_b(const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell, const Point<dim - 1> &eta);

  Tensor<1, dim>
  evaluate_Vk_integrals(const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell, const Point<dim - 1> &eta);

  Tensor<1, dim>
  evaluate_Wk_integrals(const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell, const Point<dim - 1> &eta);

  std::vector<Tensor<1, dim>>
  evaluate_VkNj_integrals(const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell, const Point<dim - 1> &eta);

  std::vector<Tensor<1, dim>>
  evaluate_WkNj_integrals(const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell, const Point<dim - 1> &eta);

  bool printstuff = false;

private:
  double                       rho_quadrature_order;
  double                       theta_quadrature_order;
  FiniteElement<dim - 1, dim> &fe;
  Mapping<dim - 1, dim> &      mapping;
  std::vector<Point<3>>        ref_vertices          = {Point<3>(0.0, 0.0, 0.0), Point<3>(1.0, 0.0, 0.0), Point<3>(0.0, 1.0, 0.0), Point<3>(1.0, 1.0, 0.0)};
  int                          ref_edge_to_vtx[4][2] = {{2, 0}, {1, 3}, {0, 1}, {3, 2}}; // Counter clockwise orientation of edges
};

#endif
