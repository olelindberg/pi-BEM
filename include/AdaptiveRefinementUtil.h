#ifndef ADAPTIVE_REFINEMENT_UTIL_H
#define ADAPTIVE_REFINEMENT_UTIL_H

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/lac/trilinos_vector.h>


#include <iostream>
#include <vector>

class AdaptiveRefinementUtil
{
public:
  typedef typename dealii::DoFHandler<2, 3>::active_cell_iterator cell_it;

  static void
  scalarFieldRanges(int                             numCellDOF,
                    const dealii::DoFHandler<2, 3> &dh,
                    const dealii::Vector<double> &  scalarField,
                    dealii::Vector<double> &        ranges)
  {
    std::vector<dealii::types::global_dof_index> indices(numCellDOF);
    for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
    {
      cell->get_dof_indices(indices);

      double minval = std::numeric_limits<double>::max();
      double maxval = -std::numeric_limits<double>::max();
      for (int j = 0; j < numCellDOF; ++j)
      {
        auto val = scalarField[indices[j]];
        minval   = std::min(val, minval);
        maxval   = std::max(val, maxval);
      } // for j in cell dofs
      ranges[cell->active_cell_index()] = maxval - minval;
    }
  }

  static void
  normalizeRanges(dealii::Vector<double> &ranges)
  {
    ranges.add(-ranges.mean_value());
    ranges /= (ranges.l2_norm() * std::sqrt(1.0 / ranges.size()));
  }

  static void
  assignRefinement(double                        errorEstimatorMax,
                   double                        aspectRatioMax,
                   const dealii::Vector<double> &error_estimator,
                   dealii::DoFHandler<2, 3> &    dh)
  {
    for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
    {
      if (error_estimator[cell->active_cell_index()] > errorEstimatorMax)
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
    } // for cell in active cells
  }
};

#endif // ADAPTIVE_REFINEMENT_UTIL_H