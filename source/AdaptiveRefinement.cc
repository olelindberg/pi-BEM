#include "../include/AdaptiveRefinement.h"
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
AdaptiveRefinement::_assignRefinement(const dealii::Vector<double> &error_estimator,
                                      dealii::DoFHandler<2, 3> &    dh)
{
  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
    if (error_estimator[cell->active_cell_index()] > _errorEstimatorMax)
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
}

void
AdaptiveRefinement::refine(unsigned int                                 np,
                           unsigned int                                 pid,
                           const dealii::Mapping<2,3>& mapping,
                           const std::vector<dealii::types::global_dof_index>& original_to_sub_wise,
                           const std::vector<dealii::types::global_dof_index>& sub_wise_to_original,
                           const dealii::FiniteElement<2, 3> &          fe,
                           const dealii::FiniteElement<2, 3> &          gradient_fe,
                           dealii::DoFHandler<2, 3> &                   dh,
                           const dealii::DoFHandler<2, 3> &             gradient_dh,
                           const dealii::TrilinosWrappers::MPI::Vector &error_vector,
                           const dealii::TrilinosWrappers::MPI::Vector &vector_gradients_solution,
                           dealii::Triangulation<2, 3> &                tria)
{





  //---------------------------------------------------------------------------
  // Adaptation to velocity potential:
  //---------------------------------------------------------------------------
    const dealii::Vector<double> error_vector_local(error_vector);
  _error_estimator_potential.reinit(tria.n_active_cells());

  std::vector<dealii::types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
      cell->get_dof_indices(local_dof_indices);

      double minval = std::numeric_limits<double>::max();
      double maxval = -std::numeric_limits<double>::max();
      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
      {
        auto val = error_vector_local[local_dof_indices[j]];
        minval   = std::min(val, minval);
        maxval   = std::max(val, maxval);
      } // for j in cell dofs

      _error_estimator_potential[cell->active_cell_index()] = maxval - minval;
      std::cout << pid << " cell " << cell->active_cell_index() << " nodes";
      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
          std::cout << " " << local_dof_indices[j];
      std::cout << "\n";
  }   // for cell in active cells

  //---------------------------------------------------------------------------
  // Substract mean and divide by standard deviation:
  //---------------------------------------------------------------------------
  _error_estimator_potential.add(-_error_estimator_potential.mean_value());
  _error_estimator_potential /=
    (_error_estimator_potential.l2_norm() * std::sqrt(1.0 / _error_estimator_potential.size()));

  //---------------------------------------------------------------------------
  // Adaptation to velocity gradient magnitude:
  //---------------------------------------------------------------------------
  dealii::TrilinosWrappers::MPI::Vector error_estimator_velocity(error_vector);
  error_estimator_velocity = 0.0;
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

  {
    auto        np_str  = std::to_string(np);
    auto        pid_str = std::to_string(pid);
    std::string filename =
      std::string("cellcenters_np_").append(np_str).append("_pid_").append(pid_str).append(".csv");
    std::cout << filename << std::endl;
    std::ofstream file(filename);
    if (file.is_open())
    {
      for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
      {
        if (cell->subdomain_id() == pid)
        {
          auto cc = cell->center();
          file << cell->active_cell_index() << " " << cc[0] << " " << cc[1] << " " << cc[2]
              << std::endl;
        } // if this cpu
      }   // for cell in active cells
    }
    else
      std::cout << "Unable to open " << filename << std::endl;
  }
  {
    std::vector<dealii::Point<3>> support_points(dh.n_dofs());
    dealii::DoFTools::map_dofs_to_support_points<2,3>(mapping, dh, support_points);

    auto        np_str  = std::to_string(np);
    auto        pid_str = std::to_string(pid);
    std::string filename =
      std::string("nodes_np_").append(np_str).append("_pid_").append(pid_str).append(".csv");
    std::cout << filename << std::endl;
    std::ofstream file(filename);
    if (file.is_open())
    {
      int cnt = 0;
      for (auto& pnt : support_points)
      {
          //file << sub_wise_to_original[cnt] << " " << pnt[0] << " " << pnt[1] << " " << pnt[2] << std::endl;
          file << cnt << " " << pnt[0] << " " << pnt[1] << " " << pnt[2] << std::endl;
        ++cnt;
      }
    }
    else
      std::cout << "Unable to open " << filename << std::endl;

  }


  //---------------------------------------------------------------------------
  // Substract mean and divide by standard deviation:
  //---------------------------------------------------------------------------
  error_estimator_velocity.add(-error_estimator_velocity.mean_value());
  error_estimator_velocity /=
    (error_estimator_velocity.l2_norm() * std::sqrt(1.0 / error_estimator_velocity.size()));

  const dealii::Vector<double> error_estimator_velocity_local(error_estimator_velocity);

  _assignRefinement(_error_estimator_potential, dh);
  //  _assignRefinement(error_estimator_velocity_local, dh);

  //  tria.prepare_coarsening_and_refinement();
  //  tria.execute_coarsening_and_refinement();
}
