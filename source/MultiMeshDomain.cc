#include "../include/MultiMeshDomain.h"
#include "../include/GridRefinementCreator.h"
#include "../include/KCS_LimfjordSetup.h"
#include "../include/ManifoldCreator.h"
#include "../include/MeshReader.h"
#include "../include/ShapesReader.h"

#include <deal.II/base/mpi.h>
#include <deal.II/grid/grid_tools.h>

#include <filesystem>


template <int dim>
MultiMeshDomain<dim>::MultiMeshDomain(MPI_Comm comm)
  : mpi_communicator(comm)
  , n_mpi_processes(dealii::Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(dealii::Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout)
{
  pcout.set_condition(this_mpi_process == 0);
}

template <int dim>
MultiMeshDomain<dim>::~MultiMeshDomain()
{}

template <int dim>
void MultiMeshDomain<dim>::read_domain(std::string input_path)
{
  std::filesystem::path root = input_path;

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
}

template <int dim>
void MultiMeshDomain<dim>::refine_and_resize(std::string input_path)
{
  pcout << "Refining and resizing ... " << std::endl;
  auto filename = std::filesystem::path(input_path).append("refinement.json").string();
  auto gridrefinement =
    GridRefinementCreator::create(filename, pcout, cad_to_projectors_tolerance_ratio * _max_tol);

  //  Writer writer;
  //  writer.save("/home/ole/dev/temp/trimeshinit.vtu", tria);

  // Do the refinement:
  // int cnt = 0;
  for (const auto &refinement : gridrefinement)
  {
    refinement->refine(tria, cad_surfaces);
    //  writer.save(
    //    std::string("/home/ole/dev/temp/trimeshrefine").append(std::to_string(cnt)).append(".vtu"),
    //    tria);
    //  ++cnt;
  }
}


template <int dim>
void MultiMeshDomain<dim>::update_triangulation()
{
  dealii::GridTools::partition_triangulation(n_mpi_processes, tria);
}


template class MultiMeshDomain<3>;
