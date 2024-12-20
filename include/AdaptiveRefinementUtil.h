#ifndef ADAPTIVE_REFINEMENT_UTIL_H
#define ADAPTIVE_REFINEMENT_UTIL_H

#include "RefinementUtil.h"

#include <deal.II/base/bounding_box.h>
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

  static void scalarFieldRanges(int                             numCellDOF,
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
      ranges[cell->active_cell_index()] = (maxval - minval);
    }
  }

  static int
  assignRefinement(double aspectRatioMax, double cellSizeMin, dealii::DoFHandler<2, 3> &dh)
  {
    int cnt = 0;
    for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
    {
      //      if (cell->user_flag_set()) // Only refine cells that are not newly coarsened
      {
        if (cell->refine_flag_set())
        {
          double l1   = cell->extent_in_direction(0);
          double l2   = cell->extent_in_direction(1);
          double lmin = std::min(l1, l2);
          if (lmin >= cellSizeMin)
          {
            RefinementUtil::aspectRatioRefinement(aspectRatioMax, cell);
            ++cnt;
          }
          else
          {
            cell->clear_refine_flag();
          }
        }
      }
    } // for cell in active cells

    return cnt;
  }

  static int assignCoarsening(double                        errorEstimatorMin,
                              const dealii::Vector<double> &error_estimator,
                              dealii::DoFHandler<2, 3> &    dh)
  {
    int cnt = 0;
    for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
    {
      if (cell->user_flag_set()) // Only coarsen cells that are not newly refined
      {
        if (error_estimator[cell->active_cell_index()] < errorEstimatorMin)
        {
          cell->set_coarsen_flag();
          ++cnt;
        }
      }
    } // for cell in active cells
    return cnt;
  }
};

#endif // ADAPTIVE_REFINEMENT_UTIL_H