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

#ifndef driver_h
#define driver_h

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
//
// #include <deal.II/lac/trilinos_vector.h>
// #include <deal.II/lac/trilinos_vector.h>
// #include <deal.II/lac/trilinos_sparse_matrix.h>
// #include <deal.II/lac/trilinos_solver.h>
// #include <deal.II/lac/trilinos_precondition.h>
//
// #include <deal.II/lac/petsc_vector.h>
// #include <deal.II/lac/petsc_parallel_vector.h>
// #include <deal.II/lac/petsc_parallel_sparse_matrix.h>
// #include <deal.II/lac/petsc_solver.h>
// #include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
// #include<deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

// And here are a few C++ standard header
// files that we will need:
#include <deal2lkit/parameter_acceptor.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "bem_fma.h"
#include "bem_problem.h"
#include "boundary_conditions.h"

#include "IPhysicalDomain.h"
#include "PiBEMSettings.h"

using namespace dealii;

// using namespace TrilinosWrappers::MPI::Vector;
// using namespace TrilinosWrappers::MPI::SparseMatrix;

/**
 * This class is in charge of organising the overall BEM simulation. It has
 * interfaces with all the other classes in order to have a complete simulation.
 */
template <int dim>
class Driver : public deal2lkit::ParameterAcceptor
{
public:
  Driver();

  ~Driver();

  /// method to declare the parameters
  /// to be read from the parameters file

  virtual void declare_parameters(ParameterHandler &prm);

  /// method to parse the needed parameters
  /// from the parameters file

  virtual void parse_parameters(ParameterHandler &prm);


  void run(std::string input_path = "", std::string output_path = "");

protected:
  ConditionalOStream pcout;

  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  std::shared_ptr<IPhysicalDomain<3>>    _physical_domain;
  std::shared_ptr<BEMProblem<3>>         bem_problem;
  std::shared_ptr<BoundaryConditions<3>> boundary_conditions;

  ParameterHandler prm;

  bool          global_refinement;
  Body          body;
  PiBEMSettings pibem_setup;

  void solve();
  void post_process(const std::string &output_path, int i);
};

#endif
