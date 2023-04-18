#include "../include/AreaRefinement.h"
#include "../include/RefinementUtil.h"
#include <BRepBndLib.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>

AreaRefinement::AreaRefinement(dealii::ConditionalOStream pcout,
                               const std::vector<int>    &manifold_id_,
                               double                     aspectRatioMax,
                               int                        levels_,
                               double                     areaMax,
                               bool                       verbose)
  : GridRefinement(pcout, verbose)
  , manifold_id(manifold_id_)
  , _aspectRatioMax(aspectRatioMax)
  , levels(levels_)
  , _areaMax(areaMax)
{}


void AreaRefinement::refine(dealii::Triangulation<2, 3> &tria)
{
  std::cout << "Refining based on area ...\n";

  for (int refineId = 0; refineId < levels; ++refineId)
  {
    int cnt = 0;
    for (dealii::Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
         cell != tria.end();
         ++cell)
    {
      if (std::find(manifold_id.begin(), manifold_id.end(), (int)cell->manifold_id()) !=
          manifold_id.end())
      {
        double area = cell->measure();
        if (area > _areaMax)
        {
          if (_verbose)
          {
            std::cout << "cell id " << cell->global_active_cell_index();
            std::cout << ", area " << area;
            std::cout << std::endl;
          }
          RefinementUtil::aspectRatioRefinement(_aspectRatioMax, cell);
          ++cnt;
        }
      }
    }
    if (cnt == 0)
      break;
    tria.prepare_coarsening_and_refinement();
    tria.execute_coarsening_and_refinement();
  }
  std::cout << "Number of global active cells: " << tria.n_global_active_cells() << std::endl;
}
