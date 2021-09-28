#include "../include/DistanceRefinement.h"

#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>

void
DistanceRefinement::refine(dealii::Triangulation<2, 3> &    tria,
                           const std::vector<TopoDS_Shape> &cad_surfaces)
{
  _pcout << "Distance refinement ... \n";

  for (int refineId = 0; refineId < levels; ++refineId)
  {
    for (dealii::Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
         cell != tria.end();
         ++cell)
    {
      if (cell->manifold_id() == manifold_id)
      {
        std::vector<dealii::Point<3>> pnts;
        for (int i = 0; i < 4; ++i)
          pnts.push_back(cell->vertex(i));

        Bnd_Box boxx;
        for (const auto pnt : pnts)
          boxx.Add(gp_Pnt(pnt[0], pnt[1], pnt[2]));

        auto dist = _box.Distance(boxx);

        if (dist < distance_max)
          cell->set_refine_flag();
      }
    }
    tria.prepare_coarsening_and_refinement();
    tria.execute_coarsening_and_refinement();
  }
}
