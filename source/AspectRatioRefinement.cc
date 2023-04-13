#include "../include/AspectRatioRefinement.h"


void AspectRatioRefinement::refine(dealii::Triangulation<2, 3> &tria)
{
  std::cout << "Refining based on aspect ratio ...\n";
  unsigned int refinedCellCounter = 1;
  unsigned int iter               = 0;
  // we repeat the aspect ratio refinement cycle until no cell has been
  // flagged for refinement, or until we reach a maximum of 10 cycles
  while (refinedCellCounter > 0 && (iter < (unsigned int)itermax))
  {
    // the refined cells counter is zeroed at the start of each cycle
    refinedCellCounter = 0;
    // we loop on the all the triangulation active cells
    for (dealii::Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
         cell != tria.end();
         ++cell)
    {
      bool refineManifold = false;
      for (const auto &id : manifold_id)
        if (cell->manifold_id() == (unsigned int)id)
          refineManifold = true;

      if (refineManifold)
      {
        // the following lines determine if the cell is more elongated
        // in its 0 or 1 direction
        unsigned int max_extent_dim = 0;
        unsigned int min_extent_dim = 1;
        if (cell->extent_in_direction(0) < cell->extent_in_direction(1))
        {
          max_extent_dim = 1;
          min_extent_dim = 0;
        }
        // we compute the extent of the cell in its maximum and minimum
        // elongation direction respectively
        double min_extent = cell->extent_in_direction(min_extent_dim);
        double max_extent = cell->extent_in_direction(max_extent_dim);

        double aspect_ratio = max_extent / min_extent;
        if (aspect_ratio > aspect_ratio_max)
        {
          if (_verbose)
          {
            std::cout << "iter " << iter;
            std::cout << ", cell id " << cell->global_active_cell_index();
            std::cout << ", aspect_ratio " << aspect_ratio;
            std::cout << ", min_extent " << min_extent;
            std::cout << ", max_extent " << max_extent;
            std::cout << std::endl;
          }
          cell->set_refine_flag(dealii::RefinementCase<2>::cut_axis(max_extent_dim));
          refinedCellCounter++;
        }
      }
    }
    // the number of cells refined in this cycle is reported before
    // proceeding with the next one
    tria.prepare_coarsening_and_refinement();
    tria.execute_coarsening_and_refinement();

    //    make_edges_conformal(_withDoubleNodes);
    ++iter;
  }
  std::cout << "Number of global active cells: " << tria.n_global_active_cells() << std::endl;
}