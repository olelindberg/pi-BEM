#include "KCS_LimfjordSetup.h"

#include "ManifoldCreator.h"
#include "MeshReader.h"
#include "MeshUtil.h"
#include "ShapesReader.h"


#include "./../tests.h"
#include "VerticalMeshProjection.h"
#include "computational_domain.h"

#include <filesystem>



int main(int argc, char **argv)
{
  std::filesystem::path root = "/home/ole/dev/projects/pi-BEM/docs/data/KCS_Limfjord";

  KCS_LimfjordSetup setup;

  auto mesh      = MeshReader::read(root, setup.mesh_inputs());
  auto shapes    = ShapesReader::read(root, setup.shape_inputs());
  auto manifolds = ManifoldCreator::make(setup.shape_inputs(), shapes, *mesh);

  // --------------------------------------------------------------------------
  // Project to manifold:
  // --------------------------------------------------------------------------
  std::vector<int> initial_projection_manifold_ids;
  for (auto input : setup.shape_inputs())
    if (input.get_mesh_projection() == "directional")
      initial_projection_manifold_ids.push_back(input.get_mesh_element_id());

  for (auto &manifold_id : initial_projection_manifold_ids)
    MeshUtil::project(*mesh, manifold_id);



  mesh->refine_global(3);

  std::string   filename0 = ("meshResult.inp");
  std::ofstream logfile0(filename0.c_str());
  logfile0 << std::setprecision(16);
  GridOut grid_out0;
  grid_out0.write_ucd(*mesh, logfile0);
  // grid_out0.write_ucd(computational_domain.tria, logfile0);
}
