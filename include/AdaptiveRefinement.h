#ifndef ADAPTIVE_REFINEMENT_H
#define ADAPTIVE_REFINEMENT_H

#include "GridRefinement.h"

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/lac/trilinos_vector.h>

#include <iostream>
#include <vector>

class AdaptiveRefinement
{
public:
  typedef typename dealii::DoFHandler<2, 3>::active_cell_iterator cell_it;

  AdaptiveRefinement(dealii::ConditionalOStream pcout,
                     MPI_Comm mpi_comm,
                     double                     errorEstimatorMax,
                     double                     aspectRatioMax);

  virtual ~AdaptiveRefinement()
  {}

  void
  refine(unsigned int                                        this_mpi_process,
         const std::unique_ptr<dealii::FiniteElement<2, 3>> &fe,
         const std::unique_ptr<dealii::FiniteElement<2, 3>> &gradient_fe,
         const dealii::DoFHandler<2, 3> &                    dh,
         const dealii::DoFHandler<2, 3> &                    gradient_dh,
         const dealii::TrilinosWrappers::MPI::Vector &       error_vector,
         const dealii::TrilinosWrappers::MPI::Vector &       vector_gradients_solution,
        const std::vector<dealii::types::global_dof_index>& original_to_sub_wise,
        const std::vector<dealii::types::global_dof_index>& sub_wise_to_original,         
         dealii::Triangulation<2, 3> &                       tria);

private:
  dealii::ConditionalOStream _pcout;
  MPI_Comm _mpi_comm;
  double _errorEstimatorMax = 0.0;
  double _aspectRatioMax    = 0.0;
};

#endif // ADAPTIVE_REFINEMENT_H