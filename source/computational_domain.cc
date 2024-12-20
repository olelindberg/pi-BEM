#include "../include/computational_domain.h"
#include "../include/ErrorMessage.h"
#include "../include/GridRefinementCreator.h"
#include "../include/Writer.h"

#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>


#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GeomLProp_SLProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Builder.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>

#include <deal2lkit/utilities.h>

#include "Teuchos_TimeMonitor.hpp"
#include "my_utilities.h"

#include <filesystem>

Teuchos::RCP<Teuchos::Time> readDomainTime = Teuchos::TimeMonitor::getNewTimer("Read domain");
Teuchos::RCP<Teuchos::Time> refineAndResizeTime =
  Teuchos::TimeMonitor::getNewTimer("Refine and resize");
Teuchos::RCP<Teuchos::Time> refineGlobalTime = Teuchos::TimeMonitor::getNewTimer("Refine global");


// @sect4{ComputationalDomain::ComputationalDomain and
// ComputationalDomain::read_parameters}
// The constructor initializes the
// variuous object in much the same
// way as done in the finite element
// programs such as step-4 or
// step-6. The only new ingredient
// here is the ParsedFunction object,
// which needs, at construction time,
// the specification of the number of
// components.
//
// For the exact solution the number
// of vector components is one, and
// no action is required since one is
// the default value for a
// ParsedFunction object. The wind,
// however, requires dim components
// to be specified. Notice that when
// declaring entries in a parameter
// file for the expression of the
// Functions::ParsedFunction, we need
// to specify the number of
// components explicitly, since the
// function
// Functions::ParsedFunction::declare_parameters
// is static, and has no knowledge of
// the number of components.
template <int dim>
ComputationalDomain<dim>::ComputationalDomain(MPI_Comm comm)
  : mpi_communicator(comm)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout)
{
  // Only output on first processor.
  pcout.set_condition(this_mpi_process == 0);
}

template <int dim>
ComputationalDomain<dim>::~ComputationalDomain()
{
  for (unsigned int i = 0; i < cad_surfaces.size(); ++i)
  {
    tria.reset_manifold(1 + i);
  }

  for (unsigned int i = 0; i < cad_curves.size(); ++i)
  {
    tria.reset_manifold(11 + i);
  }

  tria.reset_manifold(0);
}

template <int dim>
void ComputationalDomain<dim>::declare_parameters(ParameterHandler &prm)
{
  if (dim == 3)
    prm.declare_entry("Input grid name", "../grids/coarse_cube_double_nodes", Patterns::Anything());
  else
    prm.declare_entry("Input grid name", "../grids/circle", Patterns::Anything());

  prm.declare_entry("Input grid format", "inp", Patterns::Anything());

  prm.declare_entry("Input path to CAD files", "", Patterns::Anything());

  prm.declare_entry("Number of cycles", "2", Patterns::Integer());

  prm.declare_entry("Max aspect ratio", "3.5", Patterns::Double());

  prm.declare_entry("Max length", "1000.0", Patterns::Double());

  prm.declare_entry("Use iges surfaces and curves", "false", Patterns::Bool());

  prm.declare_entry("Cad tolerance to projectors tolerance ratio", "100", Patterns::Double());

  prm.declare_entry("Surface curvature adaptive refinement", "false", Patterns::Bool());

  prm.declare_entry("Cells per circle", "12", Patterns::Double());

  prm.declare_entry("Maximum number of curvature adaptive refinement cycles",
                    "5",
                    Patterns::Integer());

  prm.declare_entry("Number of global refinement to be executed before local "
                    "refinement cycle",
                    "0",
                    Patterns::Integer());

  prm.enter_subsection("Boundary Conditions ID Numbers");
  {
    prm.declare_entry("Dirichlet boundary ids", "1,110,110", Patterns::List(Patterns::Integer(0)));
    prm.declare_entry("Neumann boundary ids", "0,110,110", Patterns::List(Patterns::Integer(0)));
  }
  prm.leave_subsection();

  prm.declare_entry("Use a spheroid", "false", Patterns::Bool());

  prm.declare_entry("Axis x dimension", "2.", Patterns::Double());

  prm.declare_entry("Axis y dimension", "3.", Patterns::Double());

  prm.declare_entry("Axis z dimension", "4.", Patterns::Double());
}

template <int dim>
void ComputationalDomain<dim>::parse_parameters(ParameterHandler &prm)
{
  setup.input_grid_name   = prm.get("Input grid name");
  setup.input_grid_format = prm.get("Input grid format");
  setup.input_cad_path    = prm.get("Input path to CAD files");
  //  n_cycles                     = prm.get_integer("Number of cycles");
  // max_element_aspect_ratio     = prm.get_double("Max aspect ratio");
  // max_element_length           = prm.get_double("Max length");
  setup.use_cad_surface_and_curves = prm.get_bool("Use iges surfaces and curves");
  // surface_curvature_refinement     = prm.get_bool("Surface curvature adaptive refinement");
  // cells_per_circle = prm.get_double("Cells per circle");
  // pre_global_refinements =
  // prm.get_integer("Number of global refinement to be executed before local "
  //                "refinement cycle");
  // max_curvature_ref_cycles =
  // prm.get_integer("Maximum number of curvature adaptive refinement cycles");
  setup.cad_to_projectors_tolerance_ratio =
    prm.get_double("Cad tolerance to projectors tolerance ratio");

  // spheroid_bool   = prm.get_bool("Use a spheroid");
  // spheroid_x_axis = prm.get_double("Axis x dimension");
  // spheroid_y_axis = prm.get_double("Axis y dimension");
  // spheroid_z_axis = prm.get_double("Axis z dimension");

  prm.enter_subsection("Boundary Conditions ID Numbers");
  {
    std::vector<std::string> dirichlet_string_list =
      Utilities::split_string_list(prm.get("Dirichlet boundary ids"));
    setup.dirichlet_boundary_ids.resize(dirichlet_string_list.size());
    for (unsigned int i = 0; i < dirichlet_string_list.size(); ++i)
    {
      std::istringstream reader(dirichlet_string_list[i]);
      reader >> setup.dirichlet_boundary_ids[i];
    }

    std::vector<std::string> neumann_string_list =
      Utilities::split_string_list(prm.get("Neumann boundary ids"));
    setup.neumann_boundary_ids.resize(neumann_string_list.size());
    for (unsigned int i = 0; i < neumann_string_list.size(); ++i)
    {
      std::istringstream reader(neumann_string_list[i]);
      reader >> setup.neumann_boundary_ids[i];
    }

    // dirichlet_sur_ID1 = prm.get_integer("Dirichlet Surface 1 ID");
    // dirichlet_sur_ID2 = prm.get_integer("Dirichlet Surface 2 ID");
    // dirichlet_sur_ID3 = prm.get_integer("Dirichlet Surface 3 ID");
    // neumann_sur_ID1 = prm.get_integer("Neumann Surface 1 ID");
    // neumann_sur_ID2 = prm.get_integer("Neumann Surface 2 ID");
    // neumann_sur_ID3 = prm.get_integer("Neumann Surface 3 ID");
  }

  setup.print();

  prm.leave_subsection();
}

// @sect4{ComputationalDomain::read_domain}

// A boundary element method
// triangulation is basically the
// same as a (dim-1) dimensional
// triangulation, with the difference
// that the vertices belong to a
// (dim) dimensional space.
//
// Some of the mesh formats supported
// in deal.II use by default three
// dimensional points to describe
// meshes. These are the formats
// which are compatible with the
// boundary element method
// capabilities of deal.II. In
// particular we can use either UCD
// or GMSH formats. In both cases, we
// have to be particularly careful
// with the orientation of the mesh,
// because, unlike in the standard
// finite element case, no reordering
// or compatibility check is
// performed here.  All meshes are
// considered as oriented, because
// they are embedded in a higher
// dimensional space. (See the
// documentation of the GridIn and of
// the Triangulation for further
// details on orientation of cells in
// a triangulation.) In our case, the
// normals to the mesh are external
// to both the circle in 2d or the
// sphere in 3d.
//
// The other detail that is required
// for appropriate refinement of the
// boundary element mesh, is an
// accurate description of the
// manifold that the mesh is
// approximating. We already saw this
// several times for the boundary of
// standard finite element meshes
// (for example in step-5 and
// step-6), and here the principle
// and usage is the same, except that
// the HyperBallBoundary class takes
// an additional template parameter
// that specifies the embedding space
// dimension. The function object
// still has to be static to live at
// least as long as the triangulation
// object to which it is attached.

template <int dim>
bool ComputationalDomain<dim>::read_domain(const std::string &input_path)
{
  Teuchos::TimeMonitor localTimer(*readDomainTime);
  pcout << "\nReading domain ...\n";

  std::filesystem::path file_path =
    std::filesystem::path(input_path).append(setup.input_grid_name + "." + setup.input_grid_format);

  // Convert to absolute path
  std::filesystem::path absolute_path = std::filesystem::absolute(file_path);
  std::string           filename      = absolute_path.string();

  pcout << "Reading grid file: " << filename << std::endl;
  if (std::filesystem::is_directory(filename))
  {
    pcout << ErrorMessage::message(__FILE__, __LINE__, filename + " is a directory.");
    return false;
  }

  std::ifstream in(filename);
  if (in.is_open())
  {
    GridIn<dim - 1, dim> gi;
    gi.attach_triangulation(tria);
    if (setup.input_grid_format == "vtk")
      gi.read_vtk(in);
    else if (setup.input_grid_format == "msh")
      gi.read_msh(in);
    else if (setup.input_grid_format == "inp")
      gi.read_ucd(in, true);
    else
      Assert(false, ExcNotImplemented());

    if (setup.use_cad_surface_and_curves)
      _max_tol = this->read_cad_files_and_assign_manifold_projectors(input_path);
  }
  else
  {
    pcout << ErrorMessage::message(__FILE__, __LINE__, "Opening file " + filename + " failed.");
    return false;
  }

  return true;
}

template <int dim>
bool ComputationalDomain<dim>::read_cad_files(const std::string &input_path)
{
  pcout << "Color Files" << std::endl;
  unsigned int ii    = 1;
  bool         go_on = true;
  while (go_on == true)
  {
    std::filesystem::path cad_surface_filename =
      std::filesystem::path(input_path)
        .append(setup.input_cad_path)
        .append("Color_" + Utilities::int_to_string(ii) + ".iges");
  
    std::ifstream f(cad_surface_filename.string());
    if (f.good())
    {
      pcout << ii << "-th file exists" << std::endl;
      pcout << "Reading CAD surface file: " << cad_surface_filename << std::endl;
      TopoDS_Shape surface = OpenCASCADE::read_IGES(cad_surface_filename, 1e-3);
      cad_surfaces.push_back(surface);

      Bnd_Box boundingBox;
      BRepBndLib::Add(surface, boundingBox, Standard_False);
      double Xmin, Ymin, Zmin, Xmax, Ymax, Zmax;
      boundingBox.Get(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);
      pcout << "Bounding box: min = (" << Xmin << ", " << Ymin << ", " << Zmin << ")" << std::endl;
      pcout << "              max = (" << Xmax << ", " << Ymax << ", " << Zmax << ")" << std::endl;
      pcout << "Tolerance:    tol = " << OpenCASCADE::get_shape_tolerance(surface) << std::endl;
    }
    else
      go_on = false;
    ii++;
  }

  pcout << "Edge Files" << std::endl;
  ii    = 1;
  go_on = true;
  while (go_on == true)
  {
    std::string edge_filename =
      (setup.input_cad_path + "Curve_" + Utilities::int_to_string(ii) + ".iges");
    std::string cad_curve_filename =
      std::filesystem::path(input_path).append(edge_filename).string();
    std::ifstream f(cad_curve_filename);
    if (f.good())
    {
      pcout << ii << "-th file exists" << std::endl;
      pcout << "Reading CAD curve file: " << cad_curve_filename << std::endl;
      TopoDS_Shape curve = OpenCASCADE::read_IGES(cad_curve_filename, 1e-3);
      cad_curves.push_back(curve);

      Bnd_Box boundingBox;
      BRepBndLib::Add(curve, boundingBox, Standard_False);
      double Xmin, Ymin, Zmin, Xmax, Ymax, Zmax;
      boundingBox.Get(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);
      pcout << "Bounding box: min = (" << Xmin << ", " << Ymin << ", " << Zmin << ")" << std::endl;
      pcout << "              max = (" << Xmax << ", " << Ymax << ", " << Zmax << ")" << std::endl;
      pcout << "Tolerance:    tol = " << OpenCASCADE::get_shape_tolerance(curve) << std::endl;
    }
    else
      go_on = false;
    ii++;
  }
  return true;
} // read_cad_files

// template <>
// void
// ComputationalDomain<2>::assign_manifold_projectors(double tolerance)
// {}


template <int dim>
void ComputationalDomain<dim>::assign_manifold_projectors(double tolerance)
{
  for (unsigned int i = 0; i < cad_surfaces.size(); ++i)
  {
    pcout << "Creating normal to mesh projection manifold " << i << "\n";
    normal_to_mesh_projectors.push_back(
      std::make_shared<MyNormalToMeshProjectionManifold<2, 3>>(cad_surfaces[i], tolerance));
  }

  // FUUUUUUUUUUUUUSSSSSSSSSSSSK
  //  normal_to_mesh_projectors[cad_surfaces.size() - 1]->manifold_type = 1;

  for (unsigned int i = 0; i < cad_curves.size(); ++i)
  {
    pcout << "Creating arc length projection line manifold " << i << "\n";
    line_projectors.push_back(
      std::make_shared<OpenCASCADE::ArclengthProjectionLineManifold<2, 3>>(cad_curves[i],
                                                                           tolerance));
  }



  for (unsigned int i = 0; i < cad_surfaces.size(); ++i)
  {
    tria.set_manifold(1 + i, *normal_to_mesh_projectors[i]);
  }

  for (unsigned int i = 0; i < cad_curves.size(); ++i)
  {
    tria.set_manifold(11 + i, *line_projectors[i]);
  }

  // for (unsigned int i = 1; i < 2; ++i)
  // {
  //   tfi_manifolds.push_back(std::make_shared<TransfiniteInterpolationManifold<2, 3>>());
  //   tfi_manifolds.back()->initialize(tria);
  //   tria.set_manifold(1 + i, *tfi_manifolds.back());
  // }
  //  for (auto &tfi : tfi_manifolds)
  //    tfi->initialize(tria);
}


template <int dim>
double
ComputationalDomain<dim>::read_cad_files_and_assign_manifold_projectors(std::string input_path)
{
  double max_tol = 0;

  // Read cad files, that is Color_*.iges and Curve_*.iges files:
  this->read_cad_files(input_path);

  for (unsigned int i = 0; i < cad_surfaces.size(); ++i)
    max_tol = fmax(max_tol, OpenCASCADE::get_shape_tolerance(cad_surfaces[i]));

  for (unsigned int i = 0; i < cad_curves.size(); ++i)
    max_tol = fmax(max_tol, OpenCASCADE::get_shape_tolerance(cad_curves[i]));

  const double tolerance = setup.cad_to_projectors_tolerance_ratio * max_tol;

  this->assign_manifold_projectors(tolerance);

  return max_tol;
}

#include <memory>

template <int dim>
void ComputationalDomain<dim>::refine_and_resize(const std::string &input_path, const std::string &file_name)
{
  auto filename       = std::filesystem::path(input_path).append("refinement.json").string();
  auto gridrefinement = GridRefinementCreator::create(filename, pcout);

  // Do the refinement:
  for (const auto &refinement : gridrefinement)
    refinement->refine(tria);
}


template <int dim>
void ComputationalDomain<dim>::update_triangulation()
{
  // compute_double_vertex_cache();
  // OLE 25/3/2021 make_edges_conformal();
  // tria.execute_coarsening_and_refinement ();
  GridTools::partition_triangulation(n_mpi_processes, tria);
}

template <int dim>
void ComputationalDomain<dim>::make_edges_conformal(const bool with_double_nodes,
                                                    const bool isotropic_ref_on_opposite_side)
{
  if (with_double_nodes == false)
  {
    auto cell = tria.begin_active();
    auto endc = tria.end();
    for (cell = tria.begin_active(); cell != endc; ++cell)
    {
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
      {
        if (cell->face(f)->at_boundary())
        {
          if (cell->neighbor_index(f) > -1 && cell->neighbor_is_coarser(f))
          {
            TriaIterator<CellAccessor<dim - 1, dim>> cell_neigh = cell->neighbor(f);
            cell_neigh->set_refine_flag(RefinementCase<dim - 1>::isotropic_refinement);
          }
        }
      }
    }
    tria.prepare_coarsening_and_refinement();
    tria.execute_coarsening_and_refinement();
  }
  else
  {
    pcout << "Restoring mesh conformity on edges..." << std::endl;
    pcout << "cells before : " << tria.n_active_cells() << std::endl;
    bool to_restore = true;

    while (to_restore)
    {
      compute_double_vertex_cache();
      types::global_dof_index n_vertex     = tria.n_vertices();
      auto                    all_vertices = tria.get_vertices();

      double tol = 1e-7;
      to_restore = false;

      for (types::global_dof_index i = 0; i < n_vertex; ++i)
      {
        if (vertex_on_boundary[i] == true && double_vertex_vector[i].size() == 1)
        {
          std::vector<Point<dim>> nodes(GeometryInfo<dim - 1>::vertices_per_face);
          for (types::global_dof_index kk = 0; kk < vert_to_elems[i].size();
               ++kk) // ogni faccia ha due estremi
          {
            auto cell = vert_to_elems[i][kk];
            for (unsigned int f = 0; f < GeometryInfo<dim - 1>::faces_per_cell; ++f)
            {
              if (cell->neighbor_index(f) > -1 && cell->face(f)->at_boundary())
              {
                if (all_vertices[i].distance(cell->face(f)->vertex(0)) < tol)
                  nodes[kk] = cell->face(f)->vertex(1);
                else if (all_vertices[i].distance(cell->face(f)->vertex(1)) < tol)
                  nodes[kk] = cell->face(f)->vertex(0);
              }
            }
          }

          // we can now compute the center of the parent cell face
          Point<3> parent_face_center = 0.5 * (nodes[0] + nodes[1]);
          for (auto edgecell = edge_cells.begin(); edgecell != edge_cells.end(); ++edgecell)
          {
            for (unsigned int d = 0; d < GeometryInfo<2>::faces_per_cell; ++d)
            {
              if ((*edgecell)->face(d)->at_boundary())
              {
                if (parent_face_center.distance(
                      ((*edgecell)->face(d)->vertex(0) + (*edgecell)->face(d)->vertex(1)) / 2) <
                    tol)
                {
                  if (isotropic_ref_on_opposite_side)
                  {
                    (*edgecell)->set_refine_flag();
                    to_restore = true;
                  }
                  else
                  {
                    // otherwise, use anisotropic refinement to make
                    // edge mesh conformal

                    if ((d == 0) || (d == 1))
                      (*edgecell)->set_refine_flag(RefinementCase<2>::cut_axis(1));
                    else
                      (*edgecell)->set_refine_flag(RefinementCase<2>::cut_axis(0));
                    to_restore = true;
                  }
                }
              }
            }
          }
        }
      }

      if (to_restore)
      {
        tria.prepare_coarsening_and_refinement();
        tria.execute_coarsening_and_refinement();
      }
    }
    pcout << "cells after : " << tria.n_active_cells() << std::endl;
  }
} // make_edges_conformal

template <int dim>
void ComputationalDomain<dim>::compute_double_vertex_cache()
{
  pcout << "Computing infos for the double_vertex" << std::endl;
  double toll = 1e-7;
  double_vertex_vector.clear();
  types::global_dof_index n_vertex = tria.n_vertices();
  double_vertex_vector.resize(n_vertex);
  vertex_on_boundary.resize(n_vertex);
  std::fill(vertex_on_boundary.begin(), vertex_on_boundary.end(), false);

  auto all_vertices = tria.get_vertices();

  //---------------------------------------------------------------------------
  // Find all double vertices and connections:
  //---------------------------------------------------------------------------
  int numDoubleVertices = 0;
  for (types::global_dof_index i = 0; i < n_vertex; ++i)
  {
    for (types::global_dof_index j = 0; j < n_vertex; ++j)
    {
      if (all_vertices[i].distance(all_vertices[j]) <= toll)
      {
        double_vertex_vector[i].push_back(j);
        ++numDoubleVertices;
      }
    }
  }

  auto cell = tria.begin_active();
  auto endc = tria.end();
  vert_to_elems.clear();
  edge_cells.clear();

  for (cell = tria.begin_active(); cell != endc; ++cell)
  {
    std::vector<Point<dim>> cell_vertices(GeometryInfo<dim - 1>::vertices_per_cell);

    for (unsigned int v = 0; v < GeometryInfo<dim - 1>::vertices_per_cell; ++v)
    {
      vert_to_elems[cell->vertex_index(v)].push_back(cell);
      cell_vertices[v] = cell->vertex(v);
    }

    if (cell->at_boundary())
    {
      edge_cells.insert(cell);

      for (unsigned int f = 0; f < GeometryInfo<dim - 1>::faces_per_cell; ++f)
      {
        if (cell->face(f)->at_boundary())
        {
          for (unsigned int v = 0; v < GeometryInfo<dim - 1>::vertices_per_cell; ++v)
          {
            if (cell->face(f)->vertex(0) == cell_vertices[v])
              vertex_on_boundary[cell->vertex_index(v)] = true;
            else if (cell->face(f)->vertex(1) == cell_vertices[v])
              vertex_on_boundary[cell->vertex_index(v)] = true;
          }
        }
      }
    }
  }
}

// template class ComputationalDomain<2>;
template class ComputationalDomain<3>;
