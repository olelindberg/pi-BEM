#include "../include/MultiMeshDomain.h"
#include "../include/GridRefinementCreator.h"
#include "../include/KCS_LimfjordSetup.h"
#include "../include/ManifoldCreator.h"
#include "../include/MeshReader.h"
#include "../include/MeshUtil.h"
#include "../include/ShapesReader.h"
#include "../include/Writer.h"

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
  _mesh = std::make_shared<dealii::Triangulation<2, 3>>();

  neumann_boundary_ids = {1, 2, 3, 4};
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

  MeshReader::read(root, setup.mesh_inputs(), *_mesh);
  auto shapes    = ShapesReader::read(root, setup.shape_inputs());
  auto manifolds = ManifoldCreator::make(setup.shape_inputs(), shapes, *_mesh);

  // // --------------------------------------------------------------------------
  // // Project to manifold:
  // // --------------------------------------------------------------------------
  // std::vector<int> initial_projection_manifold_ids;
  // for (auto input : setup.shape_inputs())
  //   if (input.get_mesh_projection() == "directional")
  //     initial_projection_manifold_ids.push_back(input.get_mesh_element_id());

  // for (auto &manifold_id : initial_projection_manifold_ids)
  //   MeshUtil::project(*_mesh, manifold_id);
}

template <int dim>
void MultiMeshDomain<dim>::refine_and_resize(std::string input_path)
{
  pcout << "Refining and resizing ... " << std::endl;
  auto filename       = std::filesystem::path(input_path).append("refinement.json").string();
  auto gridrefinement = GridRefinementCreator::create(filename, pcout);

  pcout << "Saving ... " << std::endl;
  Writer writer;
  writer.save("/home/ole/dev/temp/trimeshinit.vtu", *_mesh);

  pcout << "Refining ... " << std::endl;
  // Do the refinement:
  int cnt = 0;
  for (const auto &refinement : gridrefinement)
  {
    refinement->refine(*_mesh);
    writer.save(
      std::string("/home/ole/dev/temp/trimeshrefine").append(std::to_string(cnt)).append(".vtu"),
      *_mesh);
    ++cnt;
  }
}


template <int dim>
void MultiMeshDomain<dim>::update_triangulation()
{
  dealii::GridTools::partition_triangulation(n_mpi_processes, *_mesh);
}


template class MultiMeshDomain<3>;
