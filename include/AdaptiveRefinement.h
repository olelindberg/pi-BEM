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
                     MPI_Comm                   mpi_comm,
                     double                     errorEstimatorMax,
                     double                     aspectRatioMax);

  virtual ~AdaptiveRefinement()
  {}

  void
  refine(unsigned int                                 np,
         unsigned int                                 pid,
         const dealii::FiniteElement<2, 3> &          fe,
         const dealii::FiniteElement<2, 3> &          gradient_fe,
         dealii::DoFHandler<2, 3> &                   dh,
         const dealii::DoFHandler<2, 3> &             gradient_dh,
         const dealii::TrilinosWrappers::MPI::Vector &error_vector,
         const dealii::TrilinosWrappers::MPI::Vector &vector_gradients_solution,
         dealii::Triangulation<2, 3> &                tria);

  const dealii::TrilinosWrappers::MPI::Vector &
  get_error_estimator_potential()
  {
    return _error_estimator_potential;
  }

private:
  dealii::ConditionalOStream _pcout;
  MPI_Comm                   _mpi_comm;
  double                     _errorEstimatorMax = 0.0;
  double                     _aspectRatioMax    = 0.0;

  dealii::TrilinosWrappers::MPI::Vector _error_estimator_potential;

  void
  _assignRefinement(const dealii::Vector<double> &error_estimator, dealii::DoFHandler<2, 3> &dh);
};

#endif // ADAPTIVE_REFINEMENT_H