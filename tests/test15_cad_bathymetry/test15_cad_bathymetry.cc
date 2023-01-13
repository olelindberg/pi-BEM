//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Read a file in iges format, and write it out again in the same
// format.

#include <deal.II/grid/grid_tools.h>

#include "./../tests.h"
#include "VerticalMeshProjection.h"
#include "computational_domain.h"

#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <Bnd_Box.hxx>

#include <filesystem>


class GeometryInput
{
public:
  void
  set_filename(std::string val)
  {
    _filename = val;
  }

  void
  set_scale(double val)
  {
    _scale = val;
  }

  void
  set_rotate_x(double val)
  {
    _rotate_x = val;
  }

  void
  set_rotate_y(double val)
  {
    _rotate_y = val;
  }

  void
  set_rotate_z(double val)
  {
    _rotate_z = val;
  }

  void
  set_pre_translate_x(double val)
  {
    _pre_translate_x = val;
  }

  void
  set_pre_translate_y(double val)
  {
    _pre_translate_y = val;
  }

  void
  set_pre_translate_z(double val)
  {
    _pre_translate_z = val;
  }

  void
  set_post_translate_x(double val)
  {
    _post_translate_x = val;
  }

  void
  set_post_translate_y(double val)
  {
    _post_translate_y = val;
  }

  void
  set_post_translate_z(double val)
  {
    _post_translate_z = val;
  }

  std::string
  get_filename()
  {
    return _filename;
  }

  double
  get_scale()
  {
    return _scale;
  }

  double
  get_rotate_x()
  {
    return _rotate_x;
  }

  double
  get_rotate_y()
  {
    return _rotate_y;
  }

  double
  get_rotate_z()
  {
    return _rotate_z;
  }

  double
  get_pre_translate_x()
  {
    return _pre_translate_x;
  }

  double
  get_pre_translate_y()
  {
    return _pre_translate_y;
  }

  double
  get_pre_translate_z()
  {
    return _pre_translate_z;
  }

  double
  get_post_translate_x()
  {
    return _post_translate_x;
  }

  double
  get_post_translate_y()
  {
    return _post_translate_y;
  }

  double
  get_post_translate_z()
  {
    return _post_translate_z;
  }

private:
  std::string _filename = "none";

  double _pre_translate_x = 0.0;
  double _pre_translate_y = 0.0;
  double _pre_translate_z = 0.0;

  double _scale = 1.0;

  double _rotate_x = 0.0;
  double _rotate_y = 0.0;
  double _rotate_z = 0.0;

  double _post_translate_x = 0.0;
  double _post_translate_y = 0.0;
  double _post_translate_z = 0.0;
};

class MeshProjectionCreator
{
public:
  enum class MESH_PROJECTION
  {
    ARCH_LENGTH,
    DIRECTIONAL,
    NORMAL
  };

  MeshProjectionCreator()
  {
    name_to_id.insert(std::make_pair("arch_length", MESH_PROJECTION::ARCH_LENGTH));
    name_to_id.insert(std::make_pair("directional", MESH_PROJECTION::DIRECTIONAL));
    name_to_id.insert(std::make_pair("normal", MESH_PROJECTION::NORMAL));
  }

  MESH_PROJECTION
  get_id(std::string name)
  {
    return name_to_id.at(name);
  }


private:
  std::map<std::string, MESH_PROJECTION> name_to_id;
};

class ShapeInput : public GeometryInput
{
public:
  void
  set_mesh_projection(std::string val)
  {
    _mesh_projection = val;
  }

  void
  set_mesh_element_id(int val)
  {
    _mesh_element_id = val;
  }

  std::string
  get_mesh_projection()
  {
    return _mesh_projection;
  }

  int
  get_mesh_element_id()
  {
    return _mesh_element_id;
  }

private:
  std::string _mesh_projection = "none";
  int         _mesh_element_id = 0;
};

int
main(int argc, char **argv)
{
  // Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  // // initlog();
  // MPI_Comm               mpi_communicator(MPI_COMM_WORLD);
  // ComputationalDomain<3> computational_domain(mpi_communicator);
  // deal2lkit::ParameterAcceptor::initialize("test15_cad_bathymetry.prm", "used.prm");
  // computational_domain.read_domain();


  // //  if (refine == REFINE::GLOBAL)
  // computational_domain.getTria().refine_global(2);
  // //  else if (refine == REFINE::ASPECT_RATIO)
  // //    computational_domain.aspect_ratio_refinement(4);


  // std::string   filename0 = ("meshResult.inp");
  // std::ofstream logfile0(filename0.c_str());
  // GridOut       grid_out0;
  // grid_out0.write_ucd(computational_domain.getTria(), logfile0);

  // return 0;

  std::filesystem::path root = "/home/ole/dev/projects/pi-BEM/docs/data/KCS_Limfjord";


  //---------------------------------------------------------------------------
  // Mesh input:
  //---------------------------------------------------------------------------
  std::vector<GeometryInput> mesh_inputs;
  {
    GeometryInput input;
    input.set_filename("mesh1.inp");
    //    input.set_post_translate_z(-10.8);
    mesh_inputs.push_back(input);
  }

  {
    GeometryInput input;
    input.set_filename("mesh2.inp");
    //    input.set_post_translate_z(-10.8);
    mesh_inputs.push_back(input);
  }

  //---------------------------------------------------------------------------
  // CAD geometry input:
  //---------------------------------------------------------------------------
  std::vector<ShapeInput> shape_inputs;

  {
    ShapeInput input;
    input.set_filename("Color_1.iges");
    input.set_mesh_projection("normal");
    input.set_mesh_element_id(1);
    shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("Color_2.iges");
    input.set_mesh_projection("normal");
    input.set_mesh_element_id(2);
    shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("limfjord_seabed.iges");
    input.set_mesh_projection("directional");
    input.set_mesh_element_id(3);
    //    input.set_pre_translate_x(-582883.0);
    //    input.set_pre_translate_y(-6314420.0);
    input.set_pre_translate_x(-581628.0);
    input.set_pre_translate_y(-6314958.0);
    input.set_pre_translate_z(10.8);
    input.set_rotate_z(-2.5);

    shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("Curve_1.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(11);
    shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("Curve_2.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(12);
    shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("Curve_3.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(13);
    shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("Curve_4.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(14);
    shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("Curve_5.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(15);
    //  shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("Curve_6.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(16);
    //  shape_inputs.push_back(input);
  }



  //---------------------------------------------------------------------------
  // Read meshes:
  //---------------------------------------------------------------------------
  std::vector<std::shared_ptr<Triangulation<2, 3>>> meshes;
  for (auto &input : mesh_inputs)
  {
    auto filename = std::filesystem::path(root).append(input.get_filename());
    std::cout << "Reading inp mesh file: " << filename << std::endl;
    std::ifstream file;
    file.open(filename);
    auto mesh = std::make_shared<Triangulation<2, 3>>();
    if (file.is_open())
    {
      GridIn<2, 3> grid_in;
      grid_in.attach_triangulation(*mesh);
      grid_in.read_ucd(file, true);
      file.close();
    }
    else
    {
      std::cout << filename << std::endl;
      std::cout << "Error" << std::endl;
      // Error handling
    }

    dealii::GridTools::scale(input.get_scale(), *mesh);

    dealii::GridTools::rotate(input.get_rotate_x(), 0, *mesh);
    dealii::GridTools::rotate(input.get_rotate_y(), 1, *mesh);
    dealii::GridTools::rotate(input.get_rotate_z(), 2, *mesh);

    Tensor<1, 3> translate;
    translate[0] = input.get_post_translate_x();
    translate[1] = input.get_post_translate_y();
    translate[2] = input.get_post_translate_z();

    dealii::GridTools::shift(translate, *mesh);

    meshes.push_back(mesh);
  }

  dealii::Triangulation<2, 3> mesh;
  for (auto &tmp : meshes)
    dealii::GridGenerator::merge_triangulations(
      *tmp, mesh, mesh, std::numeric_limits<double>::epsilon(), true);


  // --------------------------------------------------------------------------
  // Read shapes:
  // --------------------------------------------------------------------------
  std::vector<TopoDS_Shape> shapes;
  for (auto &input : shape_inputs)
  {
    auto filename = std::filesystem::path(root).append(input.get_filename());


    std::ifstream file(filename);
    if (file.good())
    {
      std::cout << "Reading CAD shape file: " << filename << std::endl;
      auto shape = OpenCASCADE::read_IGES(filename, 1e-3);
      shapes.push_back(shape);
    }
    else
    {
      std::cout << filename << std::endl;
      std::cout << "Error" << std::endl;
      // Error handling
    }

    {
      gp_Trsf translate;
      translate.SetTranslation(gp_Vec(input.get_pre_translate_x(),
                                      input.get_pre_translate_y(),
                                      input.get_pre_translate_z()));
      BRepBuilderAPI_Transform translationTransform(translate);
      translationTransform.Perform(shapes.back(), Standard_True);
      shapes.back() = translationTransform.Shape();
    }

    gp_Pnt  origin;
    gp_Trsf scale;
    scale.SetScale(origin, input.get_scale());
    BRepBuilderAPI_Transform scalingTransform(shapes.back(), scale);
    scalingTransform.Perform(shapes.back(), Standard_True);
    shapes.back() = scalingTransform.Shape();

    {
      gp_Trsf rotation;
      rotation.SetRotation(gp_Ax1(origin, gp_Dir(1.0, 0.0, 0.0)), input.get_rotate_x());
      BRepBuilderAPI_Transform transform(shapes.back(), rotation);
      transform.Perform(shapes.back(), Standard_True);
      shapes.back() = transform.Shape();
    }
    {
      gp_Trsf rotation;
      rotation.SetRotation(gp_Ax1(origin, gp_Dir(0.0, 1.0, 0.0)), input.get_rotate_y());
      BRepBuilderAPI_Transform transform(shapes.back(), rotation);
      transform.Perform(shapes.back(), Standard_True);
      shapes.back() = transform.Shape();
    }
    {
      gp_Trsf rotation;
      rotation.SetRotation(gp_Ax1(origin, gp_Dir(0.0, 0.0, 1.0)), input.get_rotate_z());
      BRepBuilderAPI_Transform transform(shapes.back(), rotation);
      transform.Perform(shapes.back(), Standard_True);
      shapes.back() = transform.Shape();
    }

    {
      gp_Trsf translate;
      translate.SetTranslation(gp_Vec(input.get_post_translate_x(),
                                      input.get_post_translate_y(),
                                      input.get_post_translate_z()));
      BRepBuilderAPI_Transform translationTransform(translate);
      translationTransform.Perform(shapes.back(), Standard_True);
      shapes.back() = translationTransform.Shape();
    }


    Bnd_Box boundingBox;
    BRepBndLib::Add(shapes.back(), boundingBox, Standard_False);
    double Xmin, Ymin, Zmin, Xmax, Ymax, Zmax;
    boundingBox.Get(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);
    std::cout << "Bounding box: min = (" << Xmin << ", " << Ymin << ", " << Zmin << ")"
              << std::endl;
    std::cout << "              max = (" << Xmax << ", " << Ymax << ", " << Zmax << ")"
              << std::endl;
    std::cout << "Tolerance:    tol = " << OpenCASCADE::get_shape_tolerance(shapes.back())
              << std::endl;


    std::vector<std::shared_ptr<dealii::Manifold<2, 3>>> manifolds;
    if (input.get_mesh_projection() == "directional")
    {
      std::cout << "creating directional projection" << std::endl;
      auto manifold = std::make_shared<VerticalMeshProjection<2, 3>>(
        shapes.back(), 1.0 * OpenCASCADE::get_shape_tolerance(shapes.back()));
      manifolds.push_back(manifold);
    }
    else if (input.get_mesh_projection() == "normal")
    {
      std::cout << "creating normal projection" << std::endl;
      auto manifold = std::make_shared<MyNormalToMeshProjectionManifold<2, 3>>(
        shapes.back(), 1.0 * OpenCASCADE::get_shape_tolerance(shapes.back()));
      manifolds.push_back(manifold);
    }

    else if (input.get_mesh_projection() == "arch_length")
    {
      std::cout << "creating arch length projection" << std::endl;
      auto manifold = std::make_shared<OpenCASCADE::ArclengthProjectionLineManifold<2, 3>>(
        shapes.back(), 1.0 * OpenCASCADE::get_shape_tolerance(shapes.back()));
      manifolds.push_back(manifold);
    }

    mesh.set_manifold(input.get_mesh_element_id(), *manifolds.back());
  }

  // --------------------------------------------------------------------------
  // Rotate and translate:
  // --------------------------------------------------------------------------
  // double pi    = std::acos(-1.0);
  // double angle = 0.9 * pi;
  // GridTools::rotate(angle, 2, computational_domain.tria);

  // Tensor<1, 3> translation;
  // translation[0] = 582883.0;
  // translation[1] = 6314420.0;
  // translation[2] = 0.0;
  // GridTools::shift(translation, computational_domain.tria);

  // --------------------------------------------------------------------------
  // Project to manifold:
  // --------------------------------------------------------------------------
  // auto &tria = computational_domain.tria;
  // for (auto cell = mesh.begin_active(); cell != mesh.end(); ++cell)
  // {
  //   std::vector<Point<3>> vertices(GeometryInfo<2>::vertices_per_cell);

  //   for (unsigned int i = 0; i < GeometryInfo<2>::vertices_per_cell; ++i)
  //     vertices[i] = cell->vertex(i);

  //   for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
  //     cell->vertex(v) = cell->get_manifold().project_to_manifold(vertices, cell->vertex(v));
  // }

  mesh.refine_global(5);

  std::string   filename0 = ("meshResult.inp");
  std::ofstream logfile0(filename0.c_str());
  logfile0 << std::setprecision(16);
  GridOut grid_out0;
  grid_out0.write_ucd(mesh, logfile0);
  // grid_out0.write_ucd(computational_domain.tria, logfile0);
}
