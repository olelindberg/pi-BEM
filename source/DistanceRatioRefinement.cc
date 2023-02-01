#include "../include/DistanceRatioRefinement.h"
#include "../include/RefinementUtil.h"
#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>



DistanceRatioRefinement::DistanceRatioRefinement(dealii::ConditionalOStream pcout,
                                                 const Bnd_Box &            box,
                                                 double                     aspectRatioMax,
                                                 unsigned int               manifold_id_,
                                                 int                        levels_,
                                                 double                     distanceRatioMax)
  : GridRefinement(pcout)
  , _box(box)
  , _aspectRatioMax(aspectRatioMax)
  , manifold_id(manifold_id_)
  , levels(levels_)
  , _distanceRatioMax(distanceRatioMax)
{
  auto sideLength = _box.CornerMax().XYZ() - _box.CornerMin().XYZ();
  _boxLengthMax   = std::max(std::max(sideLength.X(), sideLength.Y()), sideLength.Z());
}


void DistanceRatioRefinement::refine(dealii::Triangulation<2, 3> &tria)
{
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

        const auto dist = std::max(_box.Distance(boxx), _boxLengthMax);

        double cellLength = std::max(cell->extent_in_direction(0), cell->extent_in_direction(1));
        auto   distRatio  = cellLength / dist;

        if (distRatio > _distanceRatioMax)
          RefinementUtil::aspectRatioRefinement(_aspectRatioMax, cell);
      }
    }
    tria.prepare_coarsening_and_refinement();
    tria.execute_coarsening_and_refinement();
  }
}
