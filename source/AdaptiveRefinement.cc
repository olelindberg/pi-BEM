#include "../include/AdaptiveRefinement.h"
#include <deal.II/dofs/dof_accessor.h>

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
                           const dealii::DoFHandler<2, 3> &             dh,
                           const dealii::DoFHandler<2, 3> &             gradient_dh,
                           const dealii::TrilinosWrappers::MPI::Vector &error_vector,
                           const dealii::TrilinosWrappers::MPI::Vector &vector_gradients_solution,
                           dealii::Triangulation<2, 3> &                tria)
{
  //---------------------------------------------------------------------------
  // Adaptation to velocity potential:
  //---------------------------------------------------------------------------
  dealii::TrilinosWrappers::MPI::Vector error_estimator_potential(error_vector);
  error_estimator_potential = 0.0;

  std::vector<dealii::types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
    if (cell->subdomain_id() == pid)
    {
      cell->get_dof_indices(local_dof_indices);

      double minval = std::numeric_limits<double>::max();
      double maxval = -std::numeric_limits<double>::max();
      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
      {
        auto val = error_vector[local_dof_indices[j]];
        minval   = std::min(val, minval);
        maxval   = std::max(val, maxval);
      } // for j in cell dofs
      error_estimator_potential[cell->active_cell_index()] = maxval - minval;
    } // if this cpu
  }   // for cell in active cells

  //---------------------------------------------------------------------------
  // Substract mean and divide by standard deviation:
  //---------------------------------------------------------------------------
  error_estimator_potential.add(-error_estimator_potential.mean_value());
  error_estimator_potential /=
    (error_estimator_potential.l2_norm() * std::sqrt(1.0 / error_estimator_potential.size()));

  //---------------------------------------------------------------------------
  // Adaptation to velocity gradient magnitude:
  //---------------------------------------------------------------------------
  dealii::Vector<double>                       error_estimator_velocity(tria.n_active_cells());
  std::vector<dealii::types::global_dof_index> local_dof_indices_v(gradient_fe.dofs_per_cell);
  cell_it                                      cell_v = gradient_dh.begin_active();
  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
    if (cell->subdomain_id() == pid)
    {
      cell_v->get_dof_indices(local_dof_indices_v);
      double minval = std::numeric_limits<double>::max();
      double maxval = -std::numeric_limits<double>::max();
      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
      {
        double u   = vector_gradients_solution[local_dof_indices_v[j * 3 + 0]];
        double v   = vector_gradients_solution[local_dof_indices_v[j * 3 + 1]];
        double w   = vector_gradients_solution[local_dof_indices_v[j * 3 + 2]];
        double val = std::sqrt(u * u + v * v + w * w);
        minval     = std::min(val, minval);
        maxval     = std::max(val, maxval);
      } // for j in cell dofs
      error_estimator_velocity[cell->active_cell_index()] = maxval - minval;
    } // if this cpu
    ++cell_v;
  } // for cell in active cells

  //---------------------------------------------------------------------------
  // Substract mean and divide by standard deviation:
  //---------------------------------------------------------------------------
  error_estimator_velocity.add(-error_estimator_velocity.mean_value());
  error_estimator_velocity /=
    (error_estimator_velocity.l2_norm() * std::sqrt(1.0 / error_estimator_velocity.size()));

  const dealii::Vector<double> error_estimator_potential_local(error_estimator_potential);

  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
    if (error_estimator_potential_local[cell->active_cell_index()] > _errorEstimatorMax ||
        error_estimator_velocity[cell->active_cell_index()] > _errorEstimatorMax)
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
      if (aspect_ratio > _aspectRatioMax)
        cell->set_refine_flag(dealii::RefinementCase<2>::cut_axis(max_extent_dim));
      else
        cell->set_refine_flag();
    }
  } // for cell in active cells

  MPI_Barrier(MPI_COMM_WORLD);
  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();
  MPI_Barrier(MPI_COMM_WORLD);
}
