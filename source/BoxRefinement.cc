#include "../include/BoxRefinement.h"
#include "../include/RefinementUtil.h"
#include <BRepBndLib.hxx>
#include <chrono>
#include <gp_Pnt.hxx>

BoxRefinement::BoxRefinement(dealii::ConditionalOStream pcout,
                             const Bnd_Box             &box,
                             double                     aspectRatioMax,
                             double                     area_min,
                             unsigned int               manifold_id_,
                             int                        levels_)
  : GridRefinement(pcout)
  , _box(box)
  , _aspectRatioMax(aspectRatioMax)
  , _area_min(area_min)
  , manifold_id(manifold_id_)
  , levels(levels_)
{}


void BoxRefinement::refine(dealii::Triangulation<2, 3> &tria)
{
  std::cout << "Refining inside box ...\n";

  for (int refineId = 0; refineId < levels; ++refineId)
  {
    for (dealii::Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
         cell != tria.end();
         ++cell)
    {

      //std::cout << "manifold_id= " << cell->manifold_id() << std::endl;
      if (cell->manifold_id() == manifold_id)
      {
        Bnd_Box boxx;
        for (int i = 0; i < 4; ++i)
        {
          auto pnt = cell->vertex(i);
          boxx.Add(gp_Pnt(pnt[0], pnt[1], pnt[2]));
        }

        if (!_box.IsOut(boxx))
        {
          double aspect_ratio = RefinementUtil::aspectRatio(cell);
          double area         = cell->measure();
          if (aspect_ratio > _aspectRatioMax)
            RefinementUtil::aspectRatioRefinement(_aspectRatioMax, cell);
          else if (area > _area_min)
            RefinementUtil::aspectRatioRefinement(_aspectRatioMax, cell);
        }
      }
    }

    tria.prepare_coarsening_and_refinement();
    tria.execute_coarsening_and_refinement();
  }
  std::cout << "Number of global active cells: " << tria.n_global_active_cells() << std::endl;
}
