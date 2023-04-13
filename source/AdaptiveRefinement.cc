#include "../include/AdaptiveRefinement.h"
#include "../include/AdaptiveRefinementUtil.h"
#include "../include/Writer.h"
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>
#include <limits>

AdaptiveRefinement::AdaptiveRefinement(dealii::ConditionalOStream pcout,
                                       MPI_Comm                   mpi_comm,
                                       double                     aspectRatioMax,
                                       double                     cellSizeMin,
                                       int                        iterMax,
                                       double                     top_fraction_max,
                                       int                        number_of_elements_max)
  : _pcout(pcout)
  , _mpi_comm(mpi_comm)
  , _aspectRatioMax(aspectRatioMax)
  , _cellSizeMin(cellSizeMin)
  , _iterMax(iterMax)
  , _top_fraction_max(top_fraction_max)
  , _number_of_elements_max(number_of_elements_max)
{}

bool AdaptiveRefinement::refine(unsigned int                                 np,
                                unsigned int                                 pid,
                                const dealii::FiniteElement<2, 3>           &fe,
                                const dealii::FiniteElement<2, 3>           &gradient_fe,
                                dealii::DoFHandler<2, 3>                    &dh,
                                dealii::DoFHandler<2, 3>                    &gradient_dh,
                                const dealii::TrilinosWrappers::MPI::Vector &potential,
                                const dealii::TrilinosWrappers::MPI::Vector &velocity,
                                dealii::Triangulation<2, 3>                 &tria)
{
  double cellSizeMin = std::numeric_limits<double>::max();
  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
    double lmin = std::min(cell->extent_in_direction(0), cell->extent_in_direction(1));
    if (lmin < cellSizeMin)
      cellSizeMin = lmin;
  }
  cellSizeMin = std::max(_cellSizeMin, cellSizeMin);

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

  dealii::Vector<double> potential_local(potential);
  dealii::Vector<double> velocity_magnitude_local(velocity_magnitude);

  bool isrefined = false;

  for (int iter = 0; iter < _iterMax; ++iter)
  {
    //---------------------------------------------------------------------------
    // Adaptation to velocity potential:
    //---------------------------------------------------------------------------
    _error_estimator_pot.reinit(tria.n_active_cells());
    AdaptiveRefinementUtil::scalarFieldRanges(fe.dofs_per_cell,
                                              dh,
                                              potential_local,
                                              _error_estimator_pot);

    for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
      _error_estimator_pot[cell->active_cell_index()] =
        _error_estimator_pot[cell->active_cell_index()] *
        _error_estimator_pot[cell->active_cell_index()] / std::sqrt(cell->measure());

    //---------------------------------------------------------------------------
    // Adaptation to velocity gradient magnitude:
    //---------------------------------------------------------------------------
    _error_estimator_vel.reinit(tria.n_active_cells());
    AdaptiveRefinementUtil::scalarFieldRanges(fe.dofs_per_cell,
                                              dh,
                                              velocity_magnitude_local,
                                              _error_estimator_vel);

    double top_fraction = 0.0;
    double bot_fraction = 0.0;
    if (tria.n_active_cells() < _number_of_elements_max)
    {
      int refine_max = (_number_of_elements_max - tria.n_active_cells()) / 4;
      top_fraction =
        std::min(_top_fraction_max, double(refine_max) / double(tria.n_active_cells()));
    }

    dealii::GridRefinement::refine_and_coarsen_fixed_number(tria,
                                                            _error_estimator_pot,
                                                            top_fraction,
                                                            bot_fraction);

    dealii::GridRefinement::refine_and_coarsen_fixed_number(tria,
                                                            _error_estimator_vel,
                                                            top_fraction,
                                                            bot_fraction);


    //---------------------------------------------------------------------------
    // Assing refinement to cells via the dof handler:
    //---------------------------------------------------------------------------
    AdaptiveRefinementUtil::assignRefinement(_aspectRatioMax, cellSizeMin, dh);

    tria.prepare_coarsening_and_refinement();

    int num_refine_cells  = 0;
    int num_coarsen_cells = 0;
    for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
    {
      if (cell->refine_flag_set())
        ++num_refine_cells;

      if (cell->coarsen_flag_set())
        ++num_coarsen_cells;
    }

    if (num_refine_cells == 0 && num_coarsen_cells == 0)
    {
      break;
    }

    dealii::SolutionTransfer<2, dealii::Vector<double>, 3> potentialInterp(dh);
    potentialInterp.prepare_for_coarsening_and_refinement(potential_local);

    dealii::SolutionTransfer<2, dealii::Vector<double>, 3> velocityInterp(dh);
    velocityInterp.prepare_for_coarsening_and_refinement(velocity_magnitude_local);

    tria.execute_coarsening_and_refinement();

    dealii::GridTools::partition_triangulation(np, tria);

    dh.distribute_dofs(fe);
    dealii::DoFRenumbering::component_wise(dh);
    dealii::DoFRenumbering::subdomain_wise(dh);

    dealii::Vector<double> potential_old(potential_local);
    potential_local.reinit(dh.n_dofs());
    potentialInterp.interpolate(potential_old, potential_local);

    dealii::Vector<double> vel_mag_old(velocity_magnitude_local);
    velocity_magnitude_local.reinit(dh.n_dofs());
    velocityInterp.interpolate(vel_mag_old, velocity_magnitude_local);

    isrefined = true;
  }
  std::cout << "Number of global active cells: " << tria.n_global_active_cells() << std::endl;

  return isrefined;
}
