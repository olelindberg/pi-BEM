#include "../include/MultiMeshDomain.h"
#include "../include/GridRefinementCreator.h"
#include "../include/KCS_LimfjordSetup.h"
#include "../include/ManifoldCreator.h"
#include "../include/MeshReader.h"
// #include "../include/MeshUtil.h"
#include "../include/ShapesReader.h"
#include "../include/Writer.h"

#include "../include/BVH_SurfaceSelector.h"


#include <BRep_Tool.hxx>
#include <BVH_BoxSet.hxx>
#include <GeomAPI_IntCS.hxx>
// #include <BVH_Tools.hxx>
#include <BVH_Traverse.hxx>
#include <GC_MakeLine.hxx>
#include <Geom_Line.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Lin.hxx>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
// #include <deal.II/base/mpi.h>
// #include <deal.II/grid/grid_tools.h>

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
bool MultiMeshDomain<dim>::read_domain(std::string input_path)
{
  std::filesystem::path root = input_path;

  KCS_LimfjordSetup setup;

  MeshReader::read(root, setup.mesh_inputs(), *_mesh);
  _shapes    = ShapesReader::read(root, setup.shape_inputs());
  _manifolds = ManifoldCreator::make(setup.shape_inputs(), _shapes, *_mesh);

  return true;
}

template <int dim>
void MultiMeshDomain<dim>::refine_and_resize(std::string input_path)
{
  auto filename       = std::filesystem::path(input_path).append("refinement.json").string();
  auto gridrefinement = GridRefinementCreator::create(filename, pcout);

  // Do the refinement:
  for (const auto &refinement : gridrefinement)
    refinement->refine(*_mesh);
}

template <int dim>
void MultiMeshDomain<dim>::update_triangulation()
{
  dealii::GridTools::partition_triangulation(n_mpi_processes, *_mesh);
}

template <int dim>
void MultiMeshDomain<dim>::update_domain(double positionx, double positiony, double rotationz)
{
  _shapes[2].set_position(-positionx, -positiony, 10.8);
  _shapes[2].set_rotation(0.0, 0.0, -rotationz);

  KCS_LimfjordSetup setup;
  _manifolds.clear();
  _manifolds = ManifoldCreator::make(setup.shape_inputs(), _shapes, *_mesh);

  // --------------------------------------------------------------------------
  // Find interior points on the bottom:
  // --------------------------------------------------------------------------
  std::vector bottom_interior_vertex(_mesh->n_vertices(), false);
  for (auto cell = _mesh->begin_active(); cell != _mesh->end(); ++cell)
  {
    if (cell->manifold_id() == 3)
    {
      for (unsigned int i = 0; i < 4; ++i)
        bottom_interior_vertex[cell->vertex_index(i)] = true;

      for (unsigned int i = 0; i < 4; ++i)
        if (cell->line(i)->manifold_id() != 3)
          for (unsigned int j = 0; j < 2; ++j)
            bottom_interior_vertex[cell->line(i)->vertex_index(j)] = false;
    }
  }

  //---------------------------------------------------------------------------
  // Build BVH:
  //---------------------------------------------------------------------------
  BVH_BoxSet<double, 2, TopoDS_Face> bvh_box_set;
  TopExp_Explorer                    exp;
  for (exp.Init(_shapes[2].shape, TopAbs_FACE); exp.More(); exp.Next())
  {
    TopoDS_Face face = TopoDS::Face(exp.Current());

    Bnd_Box bnd_box;
    BRepBndLib::Add(face, bnd_box);

    bvh_vec_t          point_min(bnd_box.CornerMin().X(), bnd_box.CornerMin().Y());
    bvh_vec_t          point_max(bnd_box.CornerMax().X(), bnd_box.CornerMax().Y());
    BVH_Box<double, 2> bvh_box(point_min, point_max);

    bvh_box_set.Add(face, bvh_box);
  }
  bvh_box_set.Build();
  auto bvh = bvh_box_set.BVH();

  //---------------------------------------------------------------------------
  // Interpolate CAD surface to mesh:
  //---------------------------------------------------------------------------
  for (auto vtx = _mesh->begin_active_vertex(); vtx != _mesh->end_vertex(); ++vtx)
  {
    if (bottom_interior_vertex[vtx->vertex_index(0)])
    {
      if (std::fabs(vtx->vertex(0)[0] - 1100.0) < 1.0e-1)
        std::cout << "upper bnd\n";
      BVH_SurfaceSelector bvh_surface_selector;
      bvh_surface_selector.SetBVHSet(&bvh_box_set);
      bvh_surface_selector.set_point(bvh_vec_t(vtx->vertex(0)[0], vtx->vertex(0)[1]));
      if (bvh_surface_selector.Select() > 0)
      {
        gp_Pnt point(vtx->vertex(0)[0], vtx->vertex(0)[1], 1000.0);
        gp_Dir direction(0.0, 0.0, -1.0);
        gp_Lin line(point, direction);
        Handle(Geom_Line) line2 = GC_MakeLine(point, direction);
        for (auto &surface : bvh_surface_selector.get_surfaces())
        {
          GeomAPI_IntCS geom_intersect;
          geom_intersect.Perform(line2, surface);
          if (geom_intersect.IsDone() && geom_intersect.NbPoints() > 0)
          {
            double posz = geom_intersect.Point(1).Z();
            if (posz > 10.0)
              posz = 10.0;
            vtx->vertex(0)[2] = posz;
          }
        }
      }
    }
  }

  //---------------------------------------------------------------------------
  // Assign vertical vertex position to solution:
  //---------------------------------------------------------------------------
  dealii::FE_Q<2, 3> fe(1);

  dealii::DoFHandler<2, 3> dof_handler(*_mesh);
  dof_handler.distribute_dofs(fe);

  dealii::Vector<double> solution(dof_handler.n_dofs());

  std::vector<int>  vertex_id_to_dof(_mesh->n_vertices());
  std::vector<bool> visited(_mesh->n_vertices(), false);

  auto vertices = _mesh->get_vertices();
  for (auto cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
    for (int j = 0; j < 4; ++j)
      solution[cell->vertex_dof_index(j, 0)] = vertices[cell->vertex_index(j)][2];

  //---------------------------------------------------------------------------
  // Interpolate to hanging nodes:
  //---------------------------------------------------------------------------
  dealii::AffineConstraints<double> constraints;
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.distribute(solution);

  //---------------------------------------------------------------------------
  // Assign the interpolated vertical vertex positions to vertices:
  //---------------------------------------------------------------------------
  for (auto cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
    for (int j = 0; j < 4; ++j)
      vertices[cell->vertex_index(j)][2] = solution[cell->vertex_dof_index(j, 0)];
}



template class MultiMeshDomain<3>;
