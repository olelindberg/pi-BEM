#include "../include/BoxRefinement.h"
#include "../include/RefinementUtil.h"
#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>

BoxRefinement::BoxRefinement(dealii::ConditionalOStream pcout,
                             const Bnd_Box &            box,
                             double                     aspectRatioMax,
                             double                     cellSizeMin,
                             unsigned int               manifold_id_,
                             int                        levels_)
  : GridRefinement(pcout)
  , _box(box)
  , _aspectRatioMax(aspectRatioMax)
  , _cellSizeMin(cellSizeMin)
  , manifold_id(manifold_id_)
  , levels(levels_)
{}


void BoxRefinement::refine(dealii::Triangulation<2, 3> &tria)
{
  _pcout << "Box refinement ... \n";

  for (int refineId = 0; refineId < levels; ++refineId)
  {
    for (dealii::Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
         cell != tria.end();
         ++cell)
    {
      double l1   = cell->extent_in_direction(0);
      double l2   = cell->extent_in_direction(1);
      double lmin = std::min(l1, l2);
      if (lmin > _cellSizeMin && cell->manifold_id() == manifold_id)
      {
        Bnd_Box boxx;
        for (int i = 0; i < 4; ++i)
        {
          auto pnt = cell->vertex(i);
          boxx.Add(gp_Pnt(pnt[0], pnt[1], pnt[2]));
        }

        if (!_box.IsOut(boxx))
          RefinementUtil::aspectRatioRefinement(_aspectRatioMax, cell);
      }
    }
    tria.prepare_coarsening_and_refinement();
    tria.execute_coarsening_and_refinement();
  }

  std::cout << "Number of active cells after box refinement: " << tria.n_active_cells()
            << std::endl;
}
