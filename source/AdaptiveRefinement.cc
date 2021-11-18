#include "../include/AdaptiveRefinement.h"
#include "../include/AdaptiveRefinementUtil.h"
#include "../include/RefinementUtil.h"
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

AdaptiveRefinement::AdaptiveRefinement(dealii::ConditionalOStream pcout,
                                       MPI_Comm                   mpi_comm,
                                       double                     errorEstimatorMax,
                                       double                     aspectRatioMax)
  : _pcout(pcout)
  , _mpi_comm(mpi_comm)
  , _errorEstimatorMax(errorEstimatorMax)
  , _aspectRatioMax(aspectRatioMax)
{}

void
AdaptiveRefinement::refine(unsigned int                                 pid,
                           const dealii::FiniteElement<2, 3> &          fe,
                           const dealii::FiniteElement<2, 3> &          gradient_fe,
                           dealii::DoFHandler<2, 3> &                   dh,
                           const dealii::DoFHandler<2, 3> &             gradient_dh,
                           const dealii::TrilinosWrappers::MPI::Vector &potential,
                           const dealii::TrilinosWrappers::MPI::Vector &velocity,
                           dealii::Triangulation<2, 3> &                tria)
{
  //---------------------------------------------------------------------------
  // Adaptation to velocity potential:
  //---------------------------------------------------------------------------
  const dealii::Vector<double> potential_local(potential);
  _error_estimator_potential.reinit(tria.n_active_cells());
  AdaptiveRefinementUtil::scalarFieldRanges(fe.dofs_per_cell,
                                            dh,
                                            potential_local,
                                            _error_estimator_potential);
  AdaptiveRefinementUtil::normalizeRanges(_error_estimator_potential);

  //---------------------------------------------------------------------------
  // Parallel calculation of velocity magnitude:
  //---------------------------------------------------------------------------
  dealii::TrilinosWrappers::MPI::Vector velocity_magnitude(potential);
  velocity_magnitude = 0.0;
  std::vector<dealii::types::global_dof_index> scalar_indices(fe.dofs_per_cell);
  std::vector<dealii::types::global_dof_index> vector_indices(gradient_fe.dofs_per_cell);
  cell_it                                      cell_v = gradient_dh.begin_active();
  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
    if (cell->subdomain_id() == pid)
    {
      cell_v->get_dof_indices(vector_indices);
      cell->get_dof_indices(scalar_indices);
      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
      {
        double u                              = velocity[vector_indices[j * 3 + 0]];
        double v                              = velocity[vector_indices[j * 3 + 1]];
        double w                              = velocity[vector_indices[j * 3 + 2]];
        velocity_magnitude[scalar_indices[j]] = std::sqrt(u * u + v * v + w * w);
      } // for j in cell dofs
    }   // if this cpu
    ++cell_v;
  } // for cell in active cells

  //---------------------------------------------------------------------------
  // Adaptation to velocity gradient magnitude:
  //---------------------------------------------------------------------------
  const dealii::Vector<double> velocity_magnitude_local(velocity_magnitude);
  _error_estimator_velocity.reinit(tria.n_active_cells());
  AdaptiveRefinementUtil::scalarFieldRanges(fe.dofs_per_cell,
                                            dh,
                                            velocity_magnitude_local,
                                            _error_estimator_velocity);
  AdaptiveRefinementUtil::normalizeRanges(_error_estimator_velocity);

  //---------------------------------------------------------------------------
  // Assing refinement to cells via the dof handler:
  //---------------------------------------------------------------------------
  AdaptiveRefinementUtil::assignRefinement(_errorEstimatorMax,
                                           _aspectRatioMax,
                                           _error_estimator_potential,
                                           dh);

  AdaptiveRefinementUtil::assignRefinement(_errorEstimatorMax,
                                           _aspectRatioMax,
                                           _error_estimator_velocity,
                                           dh);

  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();
}
