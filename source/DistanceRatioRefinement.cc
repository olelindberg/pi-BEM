#include "../include/DistanceRatioRefinement.h"
#include "../include/RefinementUtil.h"
#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>



DistanceRatioRefinement::DistanceRatioRefinement(dealii::ConditionalOStream pcout,
                                                 const Bnd_Box             &box,
                                                 double                     aspectRatioMax,
                                                 double                     area_min,
                                                 unsigned int               manifold_id_,
                                                 int                        levels_,
                                                 double                     distanceRatioMax,
                                                 bool                       verbose)
  : GridRefinement(pcout, verbose)
  , _box(box)
  , _aspectRatioMax(aspectRatioMax)
  , _area_min(area_min)
  , manifold_id(manifold_id_)
  , levels(levels_)
  , _distanceRatioMax(distanceRatioMax)
{
  auto sideLength = _box.CornerMax().XYZ() - _box.CornerMin().XYZ();
  _boxLengthMax   = std::max(std::max(sideLength.X(), sideLength.Y()), sideLength.Z());
}


void DistanceRatioRefinement::refine(dealii::Triangulation<2, 3> &tria)
{
  std::cout << "Refining based on distance ratio ...\n";
  if (_verbose)
  {
    _pcout << "manifold_id      " << manifold_id << std::endl;
    _pcout << "levels           " << levels << std::endl;
    _pcout << "aspectRatioMax   " << _aspectRatioMax << std::endl;
    _pcout << "distanceRatioMax " << _distanceRatioMax << std::endl;
  }


  for (int refineId = 0; refineId < levels; ++refineId)
  {
    for (dealii::Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
         cell != tria.end();
         ++cell)
    {
      if (cell->manifold_id() == manifold_id)
      {
        double area = cell->measure();
        if (area > _area_min)
        {
          std::vector<dealii::Point<3>> pnts;
          for (int i = 0; i < 4; ++i)
            pnts.push_back(cell->vertex(i));

          Bnd_Box boxx;
          for (const auto pnt : pnts)
            boxx.Add(gp_Pnt(pnt[0], pnt[1], pnt[2]));

          const auto dist = _box.Distance(boxx);
          if (dist > 0.0)
          {
            double cellLength =
              std::max(cell->extent_in_direction(0), cell->extent_in_direction(1));
            auto distRatio = cellLength / dist;

            if (distRatio > _distanceRatioMax)
            {
              if (_verbose)
              {
                _pcout << "dist       " << dist << std::endl;
                _pcout << "cellLength " << cellLength << std::endl;
                _pcout << "distRatio  " << distRatio << std::endl;
              }
              RefinementUtil::aspectRatioRefinement(_aspectRatioMax, cell);
            }
          }
        }
      }
    }
    tria.prepare_coarsening_and_refinement();
    tria.execute_coarsening_and_refinement();
  }
  std::cout << "Number of global active cells: " << tria.n_global_active_cells() << std::endl;
}
