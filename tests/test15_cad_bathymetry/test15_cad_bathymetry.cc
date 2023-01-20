
#include "KCS_LimfjordSetup.h"

#include "ManifoldCreator.h"

#include "MeshReader.h"
#include "MeshUtil.h"
#include "ShapesReader.h"
#include "VerticalMeshProjection.h"
#include "computational_domain.h"

#include <filesystem>



int main(int argc, char **argv)
{
  std::filesystem::path root = "/home/ole/dev/projects/pi-BEM/docs/data/KCS_Limfjord";

  KCS_LimfjordSetup setup;

  dealii::Triangulation<2, 3> mesh;
  MeshReader::read(root, setup.mesh_inputs(), mesh);
  auto shapes    = ShapesReader::read(root, setup.shape_inputs());
  auto manifolds = ManifoldCreator::make(setup.shape_inputs(), shapes, mesh);

  mesh.refine_global(3);

  shapes[2].set_rotation(0.0, 0.0, 0.0);
  MeshUtil::project(mesh, 3);

  std::string   filename0 = ("meshResult.inp");
  std::ofstream logfile0(filename0.c_str());
  logfile0 << std::setprecision(16);
  GridOut grid_out0;
  grid_out0.write_ucd(mesh, logfile0);
}
