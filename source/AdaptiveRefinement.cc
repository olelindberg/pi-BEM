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
                                       double                     potentialErrorEstimatorMax,
                                       double                     velocityErrorEstimatorMax,
                                       double                     aspectRatioMax,
                                       double                     cellSizeMin,
                                       int                        iterMax)
  : _pcout(pcout)
  , _mpi_comm(mpi_comm)
  , _potentialErrorEstimatorMax(potentialErrorEstimatorMax)
  , _velocityErrorEstimatorMax(velocityErrorEstimatorMax)
  , _aspectRatioMax(aspectRatioMax)
  , _cellSizeMin(cellSizeMin)
  , _iterMax(iterMax)
{}

bool
AdaptiveRefinement::refine(unsigned int                                 np,
                           unsigned int                                 pid,
                           const dealii::FiniteElement<2, 3> &          fe,
                           const dealii::FiniteElement<2, 3> &          gradient_fe,
                           dealii::DoFHandler<2, 3> &                   dh,
                           dealii::DoFHandler<2, 3> &                   gradient_dh,
                           const dealii::TrilinosWrappers::MPI::Vector &potential,
                           const dealii::TrilinosWrappers::MPI::Vector &velocity,
                           dealii::Triangulation<2, 3> &                tria)
{
  double cellSizeMin = std::numeric_limits<double>::max();
  for (cell_it cell = dh.begin_active(); cell != dh.end(); ++cell)
  {
    double lmin = std::min(cell->extent_in_direction(0), cell->extent_in_direction(1));
    if (lmin < cellSizeMin)
      cellSizeMin = lmin;
  }
  std::cout << "Cell size min " << cellSizeMin << std::endl;
  cellSizeMin = std::max(_cellSizeMin, cellSizeMin);
  std::cout << "Cell size min " << cellSizeMin << std::endl;


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
    _error_estimator_potential.reinit(tria.n_active_cells());
    AdaptiveRefinementUtil::scalarFieldRanges(fe.dofs_per_cell,
                                              dh,
                                              potential_local,
                                              _error_estimator_potential);
    AdaptiveRefinementUtil::normalizeRanges(_error_estimator_potential);

    //---------------------------------------------------------------------------
    // Adaptation to velocity gradient magnitude:
    //---------------------------------------------------------------------------
    _error_estimator_velocity.reinit(tria.n_active_cells());
    AdaptiveRefinementUtil::scalarFieldRanges(fe.dofs_per_cell,
                                              dh,
                                              velocity_magnitude_local,
                                              _error_estimator_velocity);
    AdaptiveRefinementUtil::normalizeRanges(_error_estimator_velocity);

    //---------------------------------------------------------------------------
    // Assing refinement to cells via the dof handler:
    //---------------------------------------------------------------------------
    int numCellsPot = AdaptiveRefinementUtil::assignRefinement(
      _potentialErrorEstimatorMax, _aspectRatioMax, cellSizeMin, _error_estimator_potential, dh);
    int numCellsVel = AdaptiveRefinementUtil::assignRefinement(
      _velocityErrorEstimatorMax, _aspectRatioMax, cellSizeMin, _error_estimator_velocity, dh);

    if (numCellsPot == 0 && numCellsVel == 0)
    {
      break;
    }

    dealii::SolutionTransfer<2, dealii::Vector<double>, dealii::DoFHandler<2, 3>> scalarInterp(dh);

    tria.prepare_coarsening_and_refinement();

    scalarInterp.prepare_for_pure_refinement();

    tria.execute_coarsening_and_refinement();

    dealii::GridTools::partition_triangulation(np, tria);

    dh.distribute_dofs(fe);
    dealii::DoFRenumbering::component_wise(dh);
    dealii::DoFRenumbering::subdomain_wise(dh);

    dealii::Vector<double> potential_old(potential_local);
    potential_local.reinit(dh.n_dofs());
    scalarInterp.refine_interpolate(potential_old, potential_local);

    dealii::Vector<double> vel_mag_old(velocity_magnitude_local);
    velocity_magnitude_local.reinit(dh.n_dofs());
    scalarInterp.refine_interpolate(vel_mag_old, velocity_magnitude_local);

    Writer      writer;
    std::string postfix;
    if (writeCount < 10)
      postfix = "00000";
    else if (writeCount < 100)
      postfix = "0000";
    else if (writeCount < 1000)
      postfix = "000";
    else if (writeCount < 10000)
      postfix = "00";
    else if (writeCount < 100000)
      postfix = "0";
    else if (writeCount < 1000000)
      postfix = "";

    writer.save(std::string("/home/ole/dev/temp/trimesh")
                  .append(postfix)
                  .append(std::to_string(writeCount))
                  .append(".vtu"),
                tria);
    ++writeCount;


    _pcout << "Number of cells refinement by potential error estimator: " << numCellsPot
           << std::endl;
    _pcout << "Number of cells refinement by velocity error estimator: " << numCellsVel
           << std::endl;
    _pcout << "Number of cells after adaptive refinement: " << tria.n_active_cells() << std::endl;

    isrefined = true;
  }

  return isrefined;
}
