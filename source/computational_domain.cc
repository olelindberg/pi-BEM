

#include "../include/computational_domain.h"

#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>

#include <boost/filesystem.hpp>

#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <deal2lkit/utilities.h>

#include "Teuchos_TimeMonitor.hpp"
#include "my_utilities.h"

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
  , used_spherical_manifold(false)
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
void
ComputationalDomain<dim>::declare_parameters(ParameterHandler &prm)
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
void
ComputationalDomain<dim>::parse_parameters(ParameterHandler &prm)
{
  input_grid_name              = prm.get("Input grid name");
  input_grid_format            = prm.get("Input grid format");
  input_cad_path               = prm.get("Input path to CAD files");
  n_cycles                     = prm.get_integer("Number of cycles");
  max_element_aspect_ratio     = prm.get_double("Max aspect ratio");
  max_element_length           = prm.get_double("Max length");
  use_cad_surface_and_curves   = prm.get_bool("Use iges surfaces and curves");
  surface_curvature_refinement = prm.get_bool("Surface curvature adaptive refinement");
  cells_per_circle             = prm.get_double("Cells per circle");
  pre_global_refinements =
    prm.get_integer("Number of global refinement to be executed before local "
                    "refinement cycle");
  max_curvature_ref_cycles =
    prm.get_integer("Maximum number of curvature adaptive refinement cycles");
  cad_to_projectors_tolerance_ratio = prm.get_double("Cad tolerance to projectors tolerance ratio");

  spheroid_bool   = prm.get_bool("Use a spheroid");
  spheroid_x_axis = prm.get_double("Axis x dimension");
  spheroid_y_axis = prm.get_double("Axis y dimension");
  spheroid_z_axis = prm.get_double("Axis z dimension");

  prm.enter_subsection("Boundary Conditions ID Numbers");
  {
    std::vector<std::string> dirichlet_string_list =
      Utilities::split_string_list(prm.get("Dirichlet boundary ids"));
    dirichlet_boundary_ids.resize(dirichlet_string_list.size());
    for (unsigned int i = 0; i < dirichlet_string_list.size(); ++i)
    {
      std::istringstream reader(dirichlet_string_list[i]);
      reader >> dirichlet_boundary_ids[i];
    }

    std::vector<std::string> neumann_string_list =
      Utilities::split_string_list(prm.get("Neumann boundary ids"));
    neumann_boundary_ids.resize(neumann_string_list.size());
    for (unsigned int i = 0; i < neumann_string_list.size(); ++i)
    {
      std::istringstream reader(neumann_string_list[i]);
      reader >> neumann_boundary_ids[i];
    }

    // dirichlet_sur_ID1 = prm.get_integer("Dirichlet Surface 1 ID");
    // dirichlet_sur_ID2 = prm.get_integer("Dirichlet Surface 2 ID");
    // dirichlet_sur_ID3 = prm.get_integer("Dirichlet Surface 3 ID");
    // neumann_sur_ID1 = prm.get_integer("Neumann Surface 1 ID");
    // neumann_sur_ID2 = prm.get_integer("Neumann Surface 2 ID");
    // neumann_sur_ID3 = prm.get_integer("Neumann Surface 3 ID");
  }
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
void
ComputationalDomain<dim>::read_domain(std::string input_path)
{
  Teuchos::TimeMonitor localTimer(*readDomainTime);
  std::cout << "Reading domain ...\n";

  std::string grid_filename =
    boost::filesystem::path(input_path).append(input_grid_name + "." + input_grid_format).string();
  std::cout << "Reading grid file: " << grid_filename << std::endl;
  std::ifstream in;
  in.open(grid_filename);
  GridIn<dim - 1, dim> gi;
  gi.attach_triangulation(tria);
  if (input_grid_format == "vtk")
    gi.read_vtk(in);
  else if (input_grid_format == "msh")
    gi.read_msh(in);
  else if (input_grid_format == "inp")
    gi.read_ucd(in, true);
  else
    Assert(false, ExcNotImplemented());

  // GridTools::copy_material_to_manifold_id(tria);

  if (input_grid_name == "../grids/coarse_sphere" ||
      input_grid_name == "../grids/coarse_sphere_double_nodes" ||
      input_grid_name == "../grids/circle")
  {
    manifold = new SphericalManifold<dim - 1, dim>;
    tria.set_all_manifold_ids(0);
    tria.set_manifold(0, *manifold);
    used_spherical_manifold = true;
  }
  std::cout << "Reading domain done\n";
}

// @sect4{ComputationalDomain::refine_and_resize}

// This function globally refines the
// mesh, distributes degrees of
// freedom, and resizes matrices and
// vectors.

// template <>
// void
// ComputationalDomain<2>::refine_and_resize (const unsigned int refinement_level, std::string
// input_path)
// {
//   pcout << "Refining and resizing mesh as required" << std::endl;
//   tria.refine_global (refinement_level);
//   pcout << "We have a tria of " << tria.n_active_cells () << " cells." << std::endl;
//   GridTools::partition_triangulation (n_mpi_processes, tria);
//   std::string   filename0 = ("meshResult.inp");
//   std::ofstream logfile0 (filename0.c_str ());
//   GridOut       grid_out0;
//   grid_out0.write_ucd (tria, logfile0);
//   pcout << "...done refining and resizing mesh" << std::endl;
// }

template <int dim>
bool
ComputationalDomain<dim>::read_cad_files(std::string input_path)
{
  pcout << "Color Files" << std::endl;
  unsigned int ii    = 1;
  bool         go_on = true;
  while (go_on == true)
  {
    std::string color_filename =
      (input_cad_path + "Color_" + Utilities::int_to_string(ii) + ".iges");
    std::string cad_surface_filename =
      boost::filesystem::path(input_path).append(color_filename).string();

    std::ifstream f(cad_surface_filename);
    if (f.good())
    {
      pcout << ii << "-th file exists" << std::endl;
      std::cout << "Reading CAD surface file: " << cad_surface_filename << std::endl;
      TopoDS_Shape surface = OpenCASCADE::read_IGES(cad_surface_filename, 1e-3);
      cad_surfaces.push_back(surface);

      Bnd_Box boundingBox;
      BRepBndLib::Add(surface, boundingBox, Standard_False);
      double Xmin, Ymin, Zmin, Xmax, Ymax, Zmax;
      boundingBox.Get(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);
      std::cout << "Bounding box: min = (" << Xmin << ", " << Ymin << ", " << Zmin << ")"
                << std::endl;
      std::cout << "              max = (" << Xmax << ", " << Ymax << ", " << Zmax << ")"
                << std::endl;
      std::cout << "Tolerance:    tol = " << OpenCASCADE::get_shape_tolerance(surface) << std::endl;
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
      (input_cad_path + "Curve_" + Utilities::int_to_string(ii) + ".iges");
    std::string cad_curve_filename =
      boost::filesystem::path(input_path).append(edge_filename).string();
    std::ifstream f(cad_curve_filename);
    if (f.good())
    {
      pcout << ii << "-th file exists" << std::endl;
      std::cout << "Reading CAD curve file: " << cad_curve_filename << std::endl;
      TopoDS_Shape curve = OpenCASCADE::read_IGES(cad_curve_filename, 1e-3);
      cad_curves.push_back(curve);

      Bnd_Box boundingBox;
      BRepBndLib::Add(curve, boundingBox, Standard_False);
      double Xmin, Ymin, Zmin, Xmax, Ymax, Zmax;
      boundingBox.Get(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);
      std::cout << "Bounding box: min = (" << Xmin << ", " << Ymin << ", " << Zmin << ")"
                << std::endl;
      std::cout << "              max = (" << Xmax << ", " << Ymax << ", " << Zmax << ")"
                << std::endl;
      std::cout << "Tolerance:    tol = " << OpenCASCADE::get_shape_tolerance(curve) << std::endl;
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
void
ComputationalDomain<dim>::assign_manifold_projectors(double tolerance)
{
  for (unsigned int i = 0; i < cad_surfaces.size(); ++i)
  {
    pcout << "Creating normal to mesh projection manifold " << i << "\n";
    //      normal_to_mesh_projectors.push_back
    //      (std::make_shared<OpenCASCADE::NormalToMeshProjectionManifold<2, 3> > (cad_surfaces[i],
    //      tolerance));
    normal_to_mesh_projectors.push_back(
      std::make_shared<MyNormalToMeshProjectionManifold<2, 3>>(cad_surfaces[i], tolerance));
  }
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
}

template <int dim>
void
ComputationalDomain<dim>::refine_and_resize(const unsigned int refinement_level,
                                            std::string        input_path)
{
  Teuchos::TimeMonitor localTimer(*refineAndResizeTime);
  pcout << "Refining and resizing ... " << std::endl;
  double max_tol = 0;
  if (use_cad_surface_and_curves)
  {
    // Read cad files, that is Color_*.iges and Curve_*.iges files:
    this->read_cad_files(input_path);

    for (unsigned int i = 0; i < cad_surfaces.size(); ++i)
    {
      max_tol = fmax(max_tol, OpenCASCADE::get_shape_tolerance(cad_surfaces[i]));
    }
    for (unsigned int i = 0; i < cad_curves.size(); ++i)
    {
      max_tol = fmax(max_tol, OpenCASCADE::get_shape_tolerance(cad_curves[i]));
    }
    const double tolerance = cad_to_projectors_tolerance_ratio * max_tol;
    pcout << "Used tolerance is: " << tolerance << std::endl;

    this->assign_manifold_projectors(tolerance);
  }

  this->aspect_ratio_refinement();

  pcout << "... done refining based on element aspect ratio\n";
  pcout << "We have a tria of " << tria.n_active_cells() << " cells." << std::endl;
  {
    std::string   filename0 = ("meshResultAspect.inp");
    std::ofstream logfile0(filename0.c_str());
    GridOut       grid_out0;
    grid_out0.write_ucd(tria, logfile0);
  }
  // the following refinement cycle is based upon the original CAD
  // geometry curvature. For this reason it can be activated not only
  // when the user requires it with the surface_curvature_refinement
  // option in the input file. Of course, this is possible only if
  // also the use_cad_surface_and_curves flag is set to thrue through
  // the input file. Only in this way in fact, CAD surfaces and curves
  // are prescribed for the triangulation refinements on some of
  // its manifold ids.
  if (use_cad_surface_and_curves && surface_curvature_refinement)
  {
    pcout << "Refining based on surface curvature ...\n";
    const double tolerance          = cad_to_projectors_tolerance_ratio * max_tol;
    unsigned int refinedCellCounter = 1;
    unsigned int cycles_counter     = 0;
    // the refinement procedure is recursively repeated until no more cells
    // are flagged for refinement, or until the user specified maximum number
    // of curvature refinement cycles is reached
    while ((refinedCellCounter) && (cycles_counter < max_curvature_ref_cycles))
    {
      std::vector<double> cell_size_all;

      // the refined cells counter is zeroed at the start of each cycle
      refinedCellCounter = 0;
      // we loop on the all the triangulation active cells
      Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
      Triangulation<2, 3>::active_cell_iterator endc = tria.end();
      for (; cell != endc; ++cell)
      {
        // In the following lines, we try to come up with an estimation
        // of the cell normal. It is obtained from the average of the
        // normal to the 4 triangles in which the cell can be split using
        // the vertices and the center. The commented lines can be used
        // for checks in case something goes wrong.

        // cout<<"center: "<<cell->center()<<endl;
        // cout<<"v0: "<<cell->vertex(0)<<endl;
        // cout<<"v1: "<<cell->vertex(1)<<endl;
        // cout<<"v2: "<<cell->vertex(2)<<endl;
        // cout<<"v3: "<<cell->vertex(3)<<endl;
        Point<3> t0 = cell->vertex(0) + (-1.0) * cell->center();
        Point<3> t1 = cell->vertex(1) + (-1.0) * cell->center();
        Point<3> t2 = cell->vertex(2) + (-1.0) * cell->center();
        Point<3> t3 = cell->vertex(3) + (-1.0) * cell->center();
        // cout<<"t0: "<<t0<<endl;
        // cout<<"t1: "<<t1<<endl;
        // cout<<"t2: "<<t2<<endl;
        // cout<<"t3: "<<t3<<endl;

        Point<3> nn0(t0(1) * t1(2) - t0(2) * t1(1),
                     t0(2) * t1(0) - t0(0) * t1(2),
                     t0(0) * t1(1) - t0(1) * t1(0));
        nn0 /= nn0.norm();
        Point<3> nn1(t1(1) * t3(2) - t1(2) * t3(1),
                     t1(2) * t3(0) - t1(0) * t3(2),
                     t1(0) * t3(1) - t1(1) * t3(0));
        nn1 /= nn1.norm();
        Point<3> nn2(t3(1) * t2(2) - t3(2) * t2(1),
                     t3(2) * t2(0) - t3(0) * t2(2),
                     t3(0) * t2(1) - t3(1) * t2(0));
        nn2 /= nn2.norm();
        Point<3> nn3(t2(1) * t0(2) - t2(2) * t0(1),
                     t2(2) * t0(0) - t2(0) * t0(2),
                     t2(0) * t0(1) - t2(1) * t0(0));
        nn3 /= nn3.norm();
        Point<3> n = (nn0 + nn1 + nn2 + nn3) / 4.0;
        n /= n.norm();
        //               cout<<cell<<endl;
        //               cout<<nn0<<endl;
        //               cout<<nn1<<endl;
        //               cout<<nn2<<endl;
        //               cout<<nn3<<endl;
        //               cout<<n<<endl;
        //               cout<<cell<<"  material id:
        //               "<<int(cell->material_id())<<endl;

        // once the cell normal has beed created, we want to use it as
        // the direction of the projection onto the CAD surface first
        // though, let's check that we are using a CAD surface for the
        // refinement of the manifold_id associated with the present cell
        double cell_size;
        if (int(cell->material_id()) - 1 < (int)cad_surfaces.size())
        {
          // if so, the cad_surface associated with the present
          // manifold_id is identified...
          TopoDS_Shape neededShape = cad_surfaces[int(cell->material_id()) - 1];
          // ...and used to set up a line intersection to project the
          // cell center on the CAD surface along the direction
          // specified by the previously computed cell normal
          //          Point<3> projection = OpenCASCADE::line_intersection (neededShape,
          //          cell->center (), n, tolerance);
          Point<3> projection = my_line_intersection<3>(neededShape, cell->center(), n, tolerance);
          //                  // in correspondence with the projected
          //                  point, we ask all the
          //                  // surface differential forms
          //                  auto tup =
          //                  OpenCASCADE::closest_point_and_differential_forms(
          //                    neededShape, projection, tolerance);
          //                  // among the differential point, we
          //                  select the maximum
          //                  // absolute curvature
          //                  double max_abs_curv =
          //                  fmax(fabs(std_cxx11::get<2>(tup)),
          //                                             fabs(std_cxx11::get<3>(tup)));

          double min_curvature = 0.0;
          double max_curvature = 0.0;
          OpenCASCADE::surface_curvature(
            neededShape, projection, tolerance, min_curvature, max_curvature);

          // among the differential point, we select the maximum
          // absolute curvature
          double max_abs_curv = fmax(fabs(min_curvature), fabs(max_curvature));

          // radius is computed from the maximum absolute curvatur
          double curvature_radius = 1.0 / fmax(max_abs_curv, tolerance);

          // the target cell size is selected so that it corresponds to
          // a cells_per_circle fraction of the circumference
          // corresponding to the minimum curvature radius
          cell_size = 2.0 * dealii::numbers::PI / cells_per_circle * curvature_radius;
        }
        else
        {
          // if the cell manifold_id is not associated to a CAD
          // surface, the target cell_size is set to and extremely high
          // value, so that the cell is never refined
          cell_size = 2 * dealii::numbers::PI / cells_per_circle / tolerance;
        }
        cell_size_all.push_back(cell_size);
        // the following line si for debug puropses and should be
        // uncommented if something is not working with the refinement

        // if the cell diameter is higher than the target cell size, the
        // refinement flag is set (unless the cell is already very small
        // ---which for us means 10xtolerance)
        if ((cell->diameter() > cell_size) && (cell->diameter() > 10 * tolerance))
        {
          cell->set_refine_flag();
          refinedCellCounter++;
        }
      }
      double min_cell_size = *std::min_element(cell_size_all.begin(), cell_size_all.end());
      double max_cell_size = *std::max_element(cell_size_all.begin(), cell_size_all.end());

      std::cout << "min_cell_size: " << min_cell_size << std::endl;
      std::cout << "max_cell_size: " << max_cell_size << std::endl;

      // the number of cells to be refined in this cycle is reported, the
      // refinement is carried out and the make_edges_conformal function is
      // called to check no edge presents non comformities
      pcout << "Curvature Based Local Refinement Cycle: " << cycles_counter << " ("
            << refinedCellCounter << ")" << std::endl;
      tria.execute_coarsening_and_refinement();
      make_edges_conformal(_withDoubleNodes);
      cycles_counter++;

      // std::string filename = ( "DTMB_II_meshResult_max_curv" +
      //                         Utilities::int_to_string(int(round(cycles_counter)))
      //                         +
      //                         ".vtk" );
      // std::ofstream logfile(filename.c_str());
      // GridOut grid_out;
      // grid_out.write_vtk(tria, logfile);
      // std::string stl_filename = ( "DTMB_II_meshResult_max_curv" +
      //                         Utilities::int_to_string(int(round(cycles_counter)))
      //                         +
      //                         ".stl" );
      // SaveSTL(tria,stl_filename);
    }
    pcout << "... done refining based on surface curvature\n";
  }
  pcout << "We have a tria of " << tria.n_active_cells() << " cells." << std::endl;

  //---------------------------------------------------------------------------
  // 1) Refine locally:
  //---------------------------------------------------------------------------

  {
    for (int refineId = 0; refineId < 1; ++refineId)
    {
      Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
      Triangulation<2, 3>::active_cell_iterator endc = tria.end();
      for (; cell != endc; ++cell)
      {
        if (cell->manifold_id() == 3)
        {
          cell->set_refine_flag();
        }
      }
      tria.execute_coarsening_and_refinement();
    }
  }

  pcout << "Post local refinement number cells: " << tria.n_active_cells() << std::endl;
  {
    std::string   filename0 = ("meshResultAspect.inp");
    std::ofstream logfile0(filename0.c_str());
    GridOut       grid_out0;
    grid_out0.write_ucd(tria, logfile0);
  }



  pcout << "Refining globally ...\n";
  {
    Teuchos::TimeMonitor localTimer(*refineGlobalTime);
    tria.refine_global(refinement_level);
  }
  pcout << "... done refining globally \n";
  pcout << "We have a tria of " << tria.n_active_cells() << " cells." << std::endl;

  GridTools::partition_triangulation(n_mpi_processes, tria);
  std::string   filename0 = ("meshResult.inp");
  std::ofstream logfile0(filename0.c_str());
  GridOut       grid_out0;
  grid_out0.write_ucd(tria, logfile0);

  pcout << "... done refining and resizing\n";
}

template <int dim>
void
ComputationalDomain<dim>::aspect_ratio_refinement(const unsigned int itermax)
{
  unsigned int refinedCellCounter = 1;
  unsigned int cycles_counter     = 0;
  // we repeat the aspect ratio refinement cycle until no cell has been
  // flagged for refinement, or until we reach a maximum of 10 cycles
  pcout << "Refining based on element aspect ratio ...\n";
  while (refinedCellCounter > 0 && (cycles_counter < itermax))
  {
    // the refined cells counter is zeroed at the start of each cycle
    refinedCellCounter = 0;
    // we loop on the all the triangulation active cells
    Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
    Triangulation<2, 3>::active_cell_iterator endc = tria.end();
    for (; cell != endc; ++cell)
    {
      // the following lines determine if the cell is more elongated
      // in its 0 or 1 direction
      unsigned int max_extent_dim = 0;
      unsigned int min_extent_dim = 1;
      if (cell->extent_in_direction(0) < cell->extent_in_direction(1))
      {
        max_extent_dim = 1;
        min_extent_dim = 0;
      }
      // we compute the extent of the cell in its maximum and minimum
      // elongation direction respectively
      double min_extent = cell->extent_in_direction(min_extent_dim);
      double max_extent = cell->extent_in_direction(max_extent_dim);
      // if the aspect ratio exceeds the prescribed maximum value, the cell
      // is refined
      if (max_extent > max_element_aspect_ratio * min_extent || max_extent > max_element_length)
      {
        cell->set_refine_flag(RefinementCase<2>::cut_axis(max_extent_dim));
        refinedCellCounter++;
      }
    }
    // the number of cells refined in this cycle is reported before
    // proceeding with the next one
    pcout << "Aspect Ratio Reduction Cycle: " << cycles_counter << " (" << refinedCellCounter << ")"
          << std::endl;
    tria.execute_coarsening_and_refinement();

    // the following commented lines are here for debug puroposes: if
    // something fails during the aspect ratio reduction cycles, they
    // should be uncommented, so that a mesh file per cycle can be
    // produced to document the evolution of the mesh through the
    // refinements. If the make_edges_conformal() function
    // is suspect in creating some error, the lines can also be
    // moved after the make_edges_conformal() function is called

    // std::string filename = ( "meshIntermediateResult_" +
    //         Utilities::int_to_string(int(round(cycles_counter))) +
    //         ".inp" );
    // std::ofstream logfile(filename.c_str());
    // GridOut grid_out;
    // grid_out.write_ucd(tria, logfile);

    make_edges_conformal(_withDoubleNodes);
    cycles_counter++;
  }
}

template <int dim>
void
ComputationalDomain<dim>::conditional_refine_and_resize(const unsigned int refinement_level)
{
  pcout << "Conditionally refining and resizing mesh as required" << std::endl;

  const Point<dim> center(0, 0, 0);
  compute_double_vertex_cache();
  make_edges_conformal(_withDoubleNodes);

  for (unsigned int step = 0; step < refinement_level; ++step)
  {
    auto cell = tria.begin_active();
    auto endc = tria.end();
    for (; cell != endc; ++cell)
    {
      for (unsigned int v = 0; v < GeometryInfo<dim - 1>::vertices_per_cell; ++v)
      {
        const double distance_from_center = center.distance(cell->vertex(v));
        if (std::fabs(distance_from_center) < 1.)
        {
          cell->set_refine_flag();
          break;
        }
      }
    }
    tria.prepare_coarsening_and_refinement();
    tria.execute_coarsening_and_refinement();
    // compute_double_vertex_cache();
    make_edges_conformal(_withDoubleNodes);
  }
  update_triangulation();
}

template <int dim>
void
ComputationalDomain<dim>::update_triangulation()
{
  std::cout << "update_triangulation01 \n";
  // compute_double_vertex_cache();
  // OLE 25/3/2021 make_edges_conformal();
  std::cout << "update_triangulation02 \n";
  // tria.execute_coarsening_and_refinement ();

  pcout << "We have a tria of " << tria.n_active_cells() << " cells." << std::endl;
  GridTools::partition_triangulation(n_mpi_processes, tria);
  // std::string   filename0 = ("meshResult.inp");
  // std::ofstream logfile0(filename0.c_str());
  // GridOut       grid_out0;
  // grid_out0.write_ucd(tria, logfile0);

  // std::ostringstream filename;
  // filename << "mesh.vtu";
  // std::ofstream output(filename.str().c_str());

  //  FE_Q<dim - 1, dim> fe_dummy(1);

  // DoFHandler<dim - 1, dim> dof_handler(tria);
  // dof_handler.distribute_dofs(fe_dummy);
  // DataOut<dim - 1, DoFHandler<dim - 1, dim>> data_out;
  // data_out.attach_dof_handler(dof_handler);
  // std::vector<unsigned int> partition_int(tria.n_active_cells());
  // GridTools::get_subdomain_association(tria, partition_int);
  // const Vector<double> partitioning(partition_int.begin(),
  // partition_int.end());
  // data_out.add_data_vector(partitioning,"partitioning",DataOut<dim - 1,
  // DoFHandler<dim - 1, dim>>::type_cell_data); data_out.build_patches();
  // data_out.write_vtu(output);

  pcout << "...done refining and resizing mesh" << std::endl;
}

// template <>
// void
// ComputationalDomain<2>::make_edges_conformal(const bool with_double_nodes,
//                                              const bool isotropic_ref_on_opposite_side)
// {
//   if (with_double_nodes || !with_double_nodes || isotropic_ref_on_opposite_side ||
//       !isotropic_ref_on_opposite_side)
//     AssertThrow(true, ExcMessage("Make  edges conformal only works in 3D"));
// }

template <int dim>
void
ComputationalDomain<dim>::make_edges_conformal(const bool with_double_nodes,
                                               const bool isotropic_ref_on_opposite_side)
{
  if (with_double_nodes == false)
  {
    std::cout << "ComputationalDomain<dim>::make_edges_conformal WITHOUT double nodes ...\n";

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
    std::cout << "ComputationalDomain<dim>::make_edges_conformal WITH double nodes ...\n";
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
    pcout << "... done restoring mesh conformity" << std::endl;
  }
} // make_edges_conformal

template <int dim>
void
ComputationalDomain<dim>::compute_double_vertex_cache()
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
  std::cout << "Found " << numDoubleVertices << " double vertices." << std::endl;

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
  pcout << "done double_vertex cache" << std::endl;
}

// template class ComputationalDomain<2>;
template class ComputationalDomain<3>;
