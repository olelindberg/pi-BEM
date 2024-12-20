#include "../include/ManifoldRefinement.h"

void ManifoldRefinement::refine(dealii::Triangulation<2, 3> &tria)
{
  for (int refineId = 0; refineId < levels; ++refineId)
  {
    for (dealii::Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
         cell != tria.end();
         ++cell)
    {
      if (cell->manifold_id() == manifoldId)
        cell->set_refine_flag();
    }
    tria.execute_coarsening_and_refinement();
  }
  std::cout << "Number of global active cells: " << tria.n_global_active_cells() << std::endl;
}