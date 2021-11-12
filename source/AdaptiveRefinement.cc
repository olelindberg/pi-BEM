#include "../include/AdaptiveRefinement.h"
#include <deal.II/dofs/dof_accessor.h>


class AdaptiveRefinementUtil
{
public:
  static void
  normalizeVector(dealii::Vector<double> &vec, const MPI_Comm& mpi_comm)
  {
    //-------------------------------------------------------------------------
    // Calculate mean:
    //-------------------------------------------------------------------------
    double mean = 0.0;
    for (auto &val : vec)
      mean += val;
    mean  = dealii::Utilities::MPI::sum(mean, mpi_comm);
    mean /= vec.size();

    //-------------------------------------------------------------------------
    // Calculate standard deviation:
    //-------------------------------------------------------------------------
    double tmp = 0.0;
    for (auto &val : vec)
    {
      const auto delta = val - mean;
      tmp += delta * delta;
    }
    tmp  = dealii::Utilities::MPI::sum(tmp, mpi_comm);
    tmp /= vec.size();
    double SD = std::sqrt(tmp);

    //-------------------------------------------------------------------------
    // Normalize by vec = (vec-mean)/SD:
    //-------------------------------------------------------------------------
    for (auto &val : vec)
    {
      val -= mean;
      val /= SD;
    }

  }
};

AdaptiveRefinement::AdaptiveRefinement( dealii::ConditionalOStream pcout,
                                        MPI_Comm mpi_comm,
                                        double                     errorEstimatorMax,
                                        double                     aspectRatioMax)
  : _pcout(pcout), _mpi_comm(mpi_comm),
  _errorEstimatorMax(errorEstimatorMax)
  , _aspectRatioMax(aspectRatioMax)
{}


void
AdaptiveRefinement::refine(unsigned int                                        pid,
                           const std::unique_ptr<dealii::FiniteElement<2, 3>> &fe,
                           const std::unique_ptr<dealii::FiniteElement<2, 3>> &gradient_fe,
                           const dealii::DoFHandler<2, 3> &                    dh,
                           const dealii::DoFHandler<2, 3> &                    gradient_dh,
                           const dealii::TrilinosWrappers::MPI::Vector &       error_vector,
                           const dealii::TrilinosWrappers::MPI::Vector &vector_gradients_solution,
                           dealii::Triangulation<2, 3> &                tria)
{
  //---------------------------------------------------------------------------
  // Adaptation to velocity potential:
  //---------------------------------------------------------------------------

  std::cout << pid << " Number of active cells " << tria.n_active_cells() << std::endl;
  dealii::Vector<double> error_estimator_potential(tria.n_active_cells());

  std::vector<dealii::types::global_dof_index> local_dof_indices(fe->dofs_per_cell);

  int cnt = 0;
  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
    if (cell->subdomain_id() == pid)
    {
      cell->get_dof_indices(local_dof_indices);


      double minval = std::numeric_limits<double>::max();
      double maxval = -std::numeric_limits<double>::max();
      for (unsigned int j = 0; j < fe->dofs_per_cell; ++j)
      {
        auto val = error_vector[local_dof_indices[j]];
        minval   = std::min(val, minval);
        maxval   = std::max(val, maxval);
      } // for j in cell dofs
      error_estimator_potential[cell->active_cell_index()] = maxval - minval;
//      std::cout << pid << " cell id " << cell->active_cell_index() << std::endl;
      ++cnt;
    } // if this cpu
  }   // for cell in active cells
//  std::cout << pid << "----------------------- Number cells visited " << cnt << std::endl;

  AdaptiveRefinementUtil::normalizeVector(error_estimator_potential,_mpi_comm);

  //---------------------------------------------------------------------------
  // Adaptation to velocity gradient magnitude:
  //---------------------------------------------------------------------------
  dealii::Vector<double> error_estimator_velocity(tria.n_active_cells());

  std::vector<dealii::types::global_dof_index> local_dof_indices_v(gradient_fe->dofs_per_cell);
  cell_it                                      cell_v = gradient_dh.begin_active();

  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
    if (cell->subdomain_id() == pid)
    {
      cell_v->get_dof_indices(local_dof_indices_v);

      double minval = std::numeric_limits<double>::max();
      double maxval = -std::numeric_limits<double>::max();
      for (unsigned int j = 0; j < fe->dofs_per_cell; ++j)
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


  AdaptiveRefinementUtil::normalizeVector(error_estimator_velocity,_mpi_comm);

  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
    if (cell->subdomain_id() == pid)
    {
      if (error_estimator_potential[cell->active_cell_index()] > _errorEstimatorMax ||
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
    } // if this cpu
  }   // for cell in active cells


  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();
}
