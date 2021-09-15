// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the pi-BEM authors.
//
// This file is part of the pi-BEM library.
//
// The pi-BEM is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License version 2.1 as published by the Free Software Foundation.
// The full text of the license can be found in the file LICENSE at
// the top level of the pi-BEM distribution.
//
// Authors: Nicola Giuliani, Andrea Mola, Luca Heltai

#ifndef surface_integral_util_h
#define surface_integral_util_h
// @sect3{Include files}

// The program starts with including a bunch
// of include files that we will use in the
// various parts of the program.


#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

//#include <deal.II/lac/petsc_vector.h>
//#include <deal.II/lac/petsc_parallel_vector.h>
//#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
//#include <deal.II/lac/petsc_solver.h>
//#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>



// And here are a few C++ standard header
// files that we will need:
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/parsed_grid_refinement.h>
#include <deal2lkit/utilities.h>
#include <mpi.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include "../include/Body.h"
#include "../include/ass_leg_function.h"
#include "../include/bem_fma.h"
#include "../include/computational_domain.h"
#include "../include/constrained_matrix.h"
#include "../include/local_expansion.h"
#include "../include/multipole_expansion.h"
#include "../include/octree_block.h"

using namespace dealii;
using namespace deal2lkit;

// using namespace TrilinosWrappers;
// using namespace TrilinosWrappers::MPI;

/**
 * - surface_integral_util. This class is the core of the BEM simulation
 *   - it receives the variables vector filled in with the proper boundary
 * condition;
 *   - it creates the codimension 1 functional space setting up the BEM;
 *   - it solves the system using a preconditioned parallel GMRES solver;
 *   - it eventually interacts with the FMM accelerator.
 */
template <int dim>
class surface_integral_util : public deal2lkit::ParameterAcceptor
{
public:
  typedef typename DoFHandler<dim - 1, dim>::active_cell_iterator cell_it;
  typedef typename DoFHandler<dim - 1, dim>::active_line_iterator line_it;

  surface_integral_util(ComputationalDomain<dim> &comp_dom, const MPI_Comm comm = MPI_COMM_WORLD);

  void
  reinit();

  virtual void
  declare_parameters(ParameterHandler &prm);

  virtual void
  parse_parameters(ParameterHandler &prm);

  double
  ssurffint(const Body &body, int dimId = 1);

  ComputationalDomain<dim> &comp_dom;


  ParsedFiniteElement<dim - 1, dim>            parsed_fe;
  ParsedFiniteElement<dim - 1, dim>            parsed_gradient_fe;
  std::unique_ptr<FiniteElement<dim - 1, dim>> fe;
  std::unique_ptr<FiniteElement<dim - 1, dim>> gradient_fe;
  DoFHandler<dim - 1, dim>                     dh;
  DoFHandler<dim - 1, dim>                     gradient_dh;

  Vector<double>                         map_vector;
  std::shared_ptr<Mapping<dim - 1, dim>> mapping = std::shared_ptr<Mapping<dim - 1, dim>>();
  unsigned int                           mapping_degree;

  std::shared_ptr<Quadrature<dim - 1>> quadrature;
  std::shared_ptr<Quadrature<dim - 2>> _quadrature1d;

  std::string mapping_type;

  IndexSet this_cpu_set;

private:
  MPI_Comm _mpi_communicator;
};

#endif
