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

// We start with including a bunch
// of include files: they might be more than
// needed, we might want to check, some time
#ifndef computational_domain_h
#define computational_domain_h

#include <IPhysicalDomain.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/smartpointer.h>
// #include <deal.II/base/std_cxx11/tuple.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

// These are the headers of the opencascade support classes and
// functions. Notice that these will contain sensible data only if you
// compiled your deal.II library with support for OpenCASCADE, i.e.,
// specifying <code>-DDEAL_II_WITH_OPENCASCADE=ON</code> and
// <code>-DOPENCASCADE_DIR=/path/to/your/opencascade/installation</code>
// when calling <code>cmake</code> during deal.II configuration.
#include <deal.II/opencascade/boundary_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <TopoDS_Shape.hxx>
#include <deal2lkit/parameter_acceptor.h>
#include <mpi.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include "../include/ass_leg_function.h"
#include "../include/local_expansion.h"
#include "../include/multipole_expansion.h"
#include "../include/my_manifold_lib.h"

using namespace dealii;
using namespace deal2lkit;

class ComputationalDomainSettings
{
public:
  /// the material ID numbers in the mesh
  /// input file, for the dirichlet_nodes
  std::vector<unsigned int> dirichlet_boundary_ids;

  /// the material ID numbers in the mesh
  /// input file, for the neumann_nodes
  std::vector<unsigned int> neumann_boundary_ids;

  // flag to assess if the software will look for cad surfaces (form files
  // Color_*.iges) and curves (from files Curve_*.iges), and use such geometries
  // to refine the grid. the program will import as many curves and surfaces as
  // there are available in the present folder, and progressively associate them
  // to the manifold IDS available in the mesh file.
  //
  bool use_cad_surface_and_curves;

  // flag to require surface refinement based on CAD surface curvature. Can only
  // be activated if previous flag is true
  //
  bool surface_curvature_refinement;


  // ratio between the tolerance found in the cad files and the
  // one to be prescribed to the projectors
  double cad_to_projectors_tolerance_ratio;

  /// Strings identifying the input grid name and format
  std::string input_grid_name;
  std::string input_grid_format;
  std::string input_cad_path;

  void print()
  {
    std::cout << "\nComputational domain settings:";
    std::cout << "\ndirichlet_boundary_ids            ";
    for (auto id : dirichlet_boundary_ids)
      std::cout << id << " ";
    std::cout << "\nneumann_boundary_ids              ";
    for (auto id : neumann_boundary_ids)
      std::cout << id << " ";
    std::cout << "\nuse_cad_surface_and_curves        " << use_cad_surface_and_curves;
    std::cout << "\nsurface_curvature_refinement      " << surface_curvature_refinement;
    std::cout << "\ncad_to_projectors_tolerance_ratio " << cad_to_projectors_tolerance_ratio;
    std::cout << "\ninput_grid_name                   " << input_grid_name;
    std::cout << "\ninput_grid_format                 " << input_grid_format;
    std::cout << "\ninput_cad_path                    " << input_cad_path;
    std::cout << "\n";
  }
};


/**
 * - ComputationalDomain. This class handles, and provides to the other classes,
 * ONLY the geometry of the problem. In particular
 *  - it handles the domain decomposition using a graph partitioning tool
 * (METIS);
 *  - it reads the domain from an external file.
 */
template <int dim>
class ComputationalDomain : public deal2lkit::ParameterAcceptor, public IPhysicalDomain<dim>
{
public:
  static std::shared_ptr<ComputationalDomain<dim>> create(MPI_Comm comm = MPI_COMM_WORLD)
  {
    return std::make_shared<ComputationalDomain<dim>>(comm);
  }

  ComputationalDomain(MPI_Comm comm = MPI_COMM_WORLD);
  ~ComputationalDomain();

  virtual const Triangulation<dim - 1, dim> &getTria() const override
  {
    return tria;
  }

  virtual Triangulation<dim - 1, dim> &getTria() override
  {
    return tria;
  }

  virtual std::vector<unsigned int> get_dirichlet_boundary_ids() override
  {
    return setup.dirichlet_boundary_ids;
  };

  virtual std::vector<unsigned int> get_neumann_boundary_ids() override
  {
    return setup.neumann_boundary_ids;
  };

  /// alternative method to read initial mesh
  /// from file
  virtual bool read_domain(const std::string& input_path = "") override;

  virtual void refine_and_resize(const std::string& input_path = "", const std::string &file_name = "refinement.json") override;

  virtual void update_triangulation() override;


  /// method to refine the imported mesh
  /// according to the level requested in
  /// the parameters file

  bool read_cad_files(const std::string& input_path = "");

  void assign_manifold_projectors(double tolerance);

  // MCJ private:
  ComputationalDomainSettings setup;

  /// method to declare the parameters
  /// to be read from the parameters file

  virtual void declare_parameters(ParameterHandler &prm);

  /// method to parse the needed parameters
  /// from the parameters file

  virtual void parse_parameters(ParameterHandler &prm);



  double read_cad_files_and_assign_manifold_projectors(std::string input_path);


  /// Here are the members of the class:
  /// they are all public, as the upper level
  /// classes (bem_problem, bem_fma,
  /// free_surface) will all need to perform
  /// operations based on the greometry (and
  /// the tools to handle it) contained in
  /// this class

  /// here are some basic classes needed by
  /// the program: a triangulation, and the
  /// FiniteElement and DoFHandler classes.
  /// A second DoF handler and FiniteElement
  /// must be created in order to compute
  /// the solution gradients, which are
  /// vectorial functions

  /// The double nodes treatment that we are
  /// using requires that the computational mesh
  /// is conformal in correspondence with the
  /// sharp edges. The next function detects
  /// if after a refinement step, non conformities
  /// appear across the edges. If any is found,
  /// it is fixed by refining the cell on the
  /// coarser side of the non conformity.
  /// The optional parameter
  /// isotropic_ref_on_opposite_side can be used to
  /// decide if such cell is refined isotropically (true)
  /// or only in one direction (false), which is the one
  /// needed to obtain conformity.
  /// As said the procedure is particularly important
  /// when double nodes are used to treat sharp edges.
  /// Yet, it is also possible to employ the procedure
  /// without double nodes. This can be specified assigning
  /// false to the value of the with_double_nodes parameter.
  /// Please note that the present function is able to make
  /// conformal all the triangulations that have at most one
  /// level of non conformity between cells that are neighboring
  /// across the edge. That is why, it should be called
  /// after every single refinement cycle that is carried out
  /// in the program execution.
  void make_edges_conformal(const bool with_double_nodes              = true,
                            const bool isotropic_ref_on_opposite_side = false);

  void compute_double_vertex_cache();
  // const unsigned int fe_degree;
  // const unsigned int mapping_degree;



  Triangulation<dim - 1, dim> tria; // \todo Move to private!!!

  /// here we are just renaming the cell
  /// iterator



  // values to be imported from the
  // parameters file:



  MPI_Comm mpi_communicator;

  unsigned int n_mpi_processes;

  unsigned int this_mpi_process;

  // to deal with conformity on edges with double nodes
  std::vector<bool>                      vertex_on_boundary;
  std::vector<std::vector<unsigned int>> double_vertex_vector;
  std::map<unsigned int, std::vector<typename Triangulation<dim - 1, dim>::active_cell_iterator>>
                                                                       vert_to_elems;
  std::set<typename Triangulation<dim - 1, dim>::active_cell_iterator> edge_cells;

  Manifold<dim - 1, dim> *manifold;


  ConditionalOStream pcout;


  /// vectors containing the CAD surfaces and curves to be
  /// (optionally) used for refinement of the triangulation
  std::vector<TopoDS_Shape> cad_surfaces;
  std::vector<TopoDS_Shape> cad_curves;

  /// vectors containing the CAD surfaces and curves projectors
  /// to be (optionally) used for refinement of the triangulation
  //  std::vector<
  //    std::shared_ptr<OpenCASCADE::NormalToMeshProjectionManifold<2, 3>>>
  //    normal_to_mesh_projectors;


  std::vector<std::shared_ptr<MyNormalToMeshProjectionManifold<2, 3>>> normal_to_mesh_projectors;
  std::vector<std::shared_ptr<OpenCASCADE::ArclengthProjectionLineManifold<2, 3>>> line_projectors;

private:
  double _max_tol         = 0.0;
  bool   _withDoubleNodes = false;
};

#endif
