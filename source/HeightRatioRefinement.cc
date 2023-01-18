#include "../include/HeightRatioRefinement.h"
#include "my_utilities.h"
#include <deal.II/grid/manifold.h>

void HeightRatioRefinement::refine(dealii::Triangulation<2, 3> &tria)
{
  _pcout << "Height ratio refinement ...\n";

  for (int iter = 0; iter < itermax; ++iter)
  {
    bool isrefining = false;
    for (dealii::Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
         cell != tria.end();
         ++cell)
    {
      std::vector<dealii::Point<3>> points;
      for (int i = 0; i < 4; ++i)
        points.push_back(cell->vertex(i));
      auto   point_new    = cell->get_manifold().project_to_manifold(points, cell->center());
      double dist_max     = (point_new - cell->center()).norm();
      double cell_extent  = std::min(cell->extent_in_direction(0), cell->extent_in_direction(1));
      double height_ratio = dist_max / cell_extent;

      if (height_ratio > height_ratio_max)
      {
        cell->set_refine_flag();
        isrefining = true;
      }
    }

    if (isrefining)
    {
      tria.prepare_coarsening_and_refinement();
      tria.execute_coarsening_and_refinement();
    }
    else
      break;
  }
}
