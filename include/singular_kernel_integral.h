#ifndef singular_kernel_integral_h
#define singular_kernel_integral_h

// The most fundamental class in the library is the Triangulation class, which
// is declared here:
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
// And this is the file in which the functions are declared that create grids:
#include <deal.II/grid/grid_generator.h>

// This file contains the description of the Lagrange interpolation finite
// element:
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>


// And this file is needed for the creation of sparsity patterns of sparse
// matrices, as shown in previous examples:
#include <deal.II/dofs/dof_tools.h>

// The next two files are needed for assembling the matrix using quadrature on
// each cell. The classes declared in them will be explained below:
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_manifold.h>
// We need the following two includes for loops over cells and/or faces:
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
// Here are some functions to generate standard grids:
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_reordering.h>

// Output of grids in various graphics formats:
#include <deal.II/grid/grid_out.h>

// This is needed for C++ output:
#include <iostream>
#include <fstream>
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

template <int dim>
class SingularKernelIntegral
{
public:

SingularKernelIntegral(typename DoFHandler<dim-1, dim>::active_cell_iterator in_cell,
                       FiniteElement<dim-1, dim> &in_fe,
                       Mapping<dim-1, dim> &in_mapping,
                       Point<dim-1> &in_eta);

Tensor<1,dim> evaluate_Vk_integrals();


std::vector<Tensor<1,dim> > evaluate_VkNj_integrals();


std::vector<Tensor<1,dim> > evaluate_WkNj_integrals();


double evaluate_integral();

private:
const typename DoFHandler<dim-1, dim>::active_cell_iterator &cell;
FiniteElement<dim-1, dim> &fe;
Mapping<dim-1, dim> &mapping;
Point<dim-1> &eta;
// to be read from input file
double rho_quadrature_order = 4;
double theta_quadrature_order = 20;

};

#endif
