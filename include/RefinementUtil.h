#ifndef REFINEMENT_UTIL_H
#define REFINEMENT_UTIL_H

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <iostream>
#include <vector>

class RefinementUtil
{
public:
  template <typename cell_type>
  static void
  aspectRatioRefinement(double aspectRatioMax, cell_type &cell)
  {
    unsigned int max_extent_dim = 0;
    unsigned int min_extent_dim = 1;
    if (cell->extent_in_direction(0) < cell->extent_in_direction(1))
    {
      max_extent_dim = 1;
      min_extent_dim = 0;
    }
    double min_extent = cell->extent_in_direction(min_extent_dim);
    double max_extent = cell->extent_in_direction(max_extent_dim);

    double aspect_ratio = max_extent / min_extent;
    if (aspect_ratio > aspectRatioMax)
      cell->set_refine_flag(dealii::RefinementCase<2>::cut_axis(max_extent_dim));
    else
      cell->set_refine_flag();
  }
};

#endif // REFINEMENT_UTIL_H