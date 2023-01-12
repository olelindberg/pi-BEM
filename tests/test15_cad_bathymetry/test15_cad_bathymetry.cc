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
#include "computational_domain.h"
#include <initializer_list>
		
class GeometryInput
{
public:

  void set_filename   (std::string val){ _filename    = val;} 
  void set_scale      (double val){ _scale       = val;}
  void set_rotate_x   (double val){ _rotate_x    = val;}
  void set_rotate_y   (double val){ _rotate_y    = val;}
  void set_rotate_z   (double val){ _rotate_z    = val;}
  void set_translate_x(double val){ _translate_x = val;}
  void set_translate_y(double val){ _translate_y = val;}
  void set_translate_z(double val){ _translate_z = val;}

  std::string get_filename   (){return _filename    ;}  
  double      get_scale      (){return _scale       ;} 
  double      get_rotate_x   (){return _rotate_x    ;} 
  double      get_rotate_y   (){return _rotate_y    ;} 
  double      get_rotate_z   (){return _rotate_z    ;} 
  double      get_translate_x(){return _translate_x ;} 
  double      get_translate_y(){return _translate_y ;} 
  double      get_translate_z(){return _translate_z ;} 

private:

  std::string _filename  = "none";
  double _scale          = 1.0;
  double _rotate_x       = 0.0;
  double _rotate_y       = 0.0;
  double _rotate_z       = 0.0;
  double _translate_x    = 0.0;
  double _translate_y    = 0.0;
  double _translate_z    = 0.0;

};

class MeshProjectionCreator
{
public:

  enum class MESH_PROJECTION 
  {
    ARCH_LENGTH, DIRECTIONAL, NORMAL
  };

  MeshProjectionCreator()
  {
    name_to_id.insert(std::make_pair("arch_length",MESH_PROJECTION::ARCH_LENGTH));
    name_to_id.insert(std::make_pair("directional",MESH_PROJECTION::DIRECTIONAL));
    name_to_id.insert(std::make_pair("normal",MESH_PROJECTION::NORMAL));
  }

  MESH_PROJECTION get_id(std::string name)
  {
    return name_to_id.at(name);
  }


private:

  std::map<std::string, MESH_PROJECTION> name_to_id;

};

class ShapeInput : public GeometryInput
{
public:

  void set_mesh_projection(std::string val){_mesh_projection = val;}
  void set_mesh_element_id(int val){_mesh_element_id = val;}

  std::string get_mesh_projection(){return _mesh_projection;}
  int get_mesh_element_id(){return _mesh_element_id;}

private:

  std::string _mesh_projection = "none";
  int _mesh_element_id = 0;

};

int
main(int argc, char **argv)
{

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  // initlog();
  MPI_Comm               mpi_communicator(MPI_COMM_WORLD);
  // ComputationalDomain<3> computational_domain(mpi_communicator);
  // deal2lkit::ParameterAcceptor::initialize("test15_cad_bathymetry.prm", "used.prm");
  // computational_domain.read_domain();

  //---------------------------------------------------------------------------
  // Mesh input:
  //---------------------------------------------------------------------------
  std::vector<GeometryInput> mesh_inputs;
  {
    GeometryInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/KCS.inp");
    input.set_translate_z(-10.8);
//    mesh_inputs.push_back(input);
  }

  {
    GeometryInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/seabed_with_walls.inp");
    input.set_scale(200.0);
    input.set_rotate_z(2.7);
    input.set_translate_x(582883.0);
    input.set_translate_y(6314420.0);

    mesh_inputs.push_back(input);
  }

  //---------------------------------------------------------------------------
  // CAD geometry input:
  //---------------------------------------------------------------------------
  std::vector<ShapeInput> shape_inputs;  
  {
    ShapeInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/limfjord_seabed.iges");
    input.set_mesh_projection("directional");
    input.set_mesh_element_id(3);
    shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/KCS_hull_stb.iges");
    input.set_mesh_projection("normal");
    input.set_mesh_element_id(1);
   // shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/KCS_hull_prt.iges");
    input.set_mesh_projection("normal");
    input.set_mesh_element_id(2);
  //  shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/waterline_stb.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(11);
  //  shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/waterline_prt.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(12);
  //  shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/centerline.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(13);
   // shape_inputs.push_back(input);
  }        

  {
    ShapeInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/cross_section_center.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(14);
  //  shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/cross_section_fwd_prt.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(15);
  //  shape_inputs.push_back(input);
  }

  {
    ShapeInput input;
    input.set_filename("/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/cross_section_fwd_stb.iges");
    input.set_mesh_projection("arch_length");
    input.set_mesh_element_id(16);
  //  shape_inputs.push_back(input);
  }        



  //---------------------------------------------------------------------------
  // Read meshes:
  //---------------------------------------------------------------------------
  std::vector<std::shared_ptr<Triangulation<2,3>>> meshes;
  for (auto& input : mesh_inputs)
  {
    std::ifstream file;
    file.open(input.get_filename());
   auto mesh = std::make_shared<Triangulation<2,3>>();
    if (file.is_open())
    {
      GridIn<2,3> grid_in;
      grid_in.attach_triangulation(*mesh);
      grid_in.read_ucd(file, true);
      file.close();
    }
    else
    {
      std::cout <<input.get_filename() << std::endl;
      std::cout << "Error" << std::endl;
      // Error handling
    }

    dealii::GridTools::scale (input.get_scale(),*mesh);      

    dealii::GridTools::rotate(input.get_rotate_x(),0,*mesh);      
    dealii::GridTools::rotate(input.get_rotate_y(),1,*mesh);      
    dealii::GridTools::rotate(input.get_rotate_z(),2,*mesh);      

    Tensor<1,3> translate;
    translate[0] = input.get_translate_x();
    translate[1] = input.get_translate_y();
    translate[2] = input.get_translate_z();

    dealii::GridTools::shift(translate,*mesh);      

    meshes.push_back(mesh);
  }

  dealii::Triangulation<2,3> mesh;
  for (auto& tmp : meshes)
      dealii::GridGenerator::merge_triangulations(*tmp,mesh,mesh,std::numeric_limits<double>::epsilon(),true);


  // --------------------------------------------------------------------------
  // Read shapes:
  // --------------------------------------------------------------------------
  std::vector<TopoDS_Shape> shapes;
  for (auto& input : shape_inputs)
  {
    std::ifstream file(input.get_filename());
    if (file.good())
    {
      std::cout << "Reading CAD shape file: " << input.get_filename() << std::endl;
      auto shape = OpenCASCADE::read_IGES(input.get_filename(), 1e-3);
      shapes.push_back(shape);
    }
    else
    {
      std::cout <<input.get_filename() << std::endl;
      std::cout << "Error" << std::endl;      
      // Error handling
    }

  std::vector<std::shared_ptr<dealii::Manifold<2,3>>> manifolds;
    if (input.get_mesh_projection()=="directional")
    {
      auto manifold = std::make_shared<MyNormalToMeshProjectionManifold<2, 3>>(shapes.back(), OpenCASCADE::get_shape_tolerance(shapes.back()));
      manifold->manifold_type = 1;
      manifolds.push_back(manifold);
    }
    else if (input.get_mesh_projection()=="normal")
    {
      auto manifold = std::make_shared<MyNormalToMeshProjectionManifold<2, 3>>(shapes.back(), OpenCASCADE::get_shape_tolerance(shapes.back()));
      manifolds.push_back(manifold);
    }

    else if (input.get_mesh_projection()=="arch_length")
    {
      auto manifold = std::make_shared<OpenCASCADE::ArclengthProjectionLineManifold<2, 3>>(shapes.back(), OpenCASCADE::get_shape_tolerance(shapes.back()));
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
   for (auto cell = mesh.begin_active(); cell != mesh.end(); ++cell)
   {
     std::vector<Point<3>> vertices(GeometryInfo<2>::vertices_per_cell);
     
      for (unsigned int i = 0; i < GeometryInfo<2>::vertices_per_cell; ++i)
        vertices[i] = cell->vertex(i);

      for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
        cell->vertex(v) = cell->get_manifold().project_to_manifold(vertices, cell->vertex(v));
   }

  mesh.refine_global(3);

  std::string   filename0 = ("meshResult.inp");
  std::ofstream logfile0(filename0.c_str());
  logfile0 << std::setprecision(16);
  GridOut grid_out0;
  grid_out0.write_ucd(mesh, logfile0);
}
