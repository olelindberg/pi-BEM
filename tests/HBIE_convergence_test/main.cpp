// #include <deal.II/base/conditional_ostream.h>
// #include <deal.II/base/convergence_table.h>
// #include <deal.II/base/parameter_acceptor.h>
// #include <deal.II/base/parsed_function.h>
// #include <deal.II/base/quadrature_lib.h>
// #include <deal.II/base/quadrature_selector.h>
// #include <deal.II/base/smartpointer.h>
// #include <deal.II/base/utilities.h>

// #include <deal.II/lac/full_matrix.h>
// #include <deal.II/lac/precondition.h>
// #include <deal.II/lac/solver_control.h>
// #include <deal.II/lac/solver_gmres.h>
// #include <deal.II/lac/sparse_matrix.h>
// #include <deal.II/lac/vector.h>

// #include <deal.II/grid/grid_generator.h>
// #include <deal.II/grid/grid_in.h>
// #include <deal.II/grid/grid_out.h>
// #include <deal.II/grid/tria.h>
// #include <deal.II/grid/tria_accessor.h>
// #include <deal.II/grid/tria_iterator.h>

// #include <deal.II/dofs/dof_accessor.h>
// #include <deal.II/dofs/dof_handler.h>
// #include <deal.II/dofs/dof_renumbering.h>
// #include <deal.II/dofs/dof_tools.h>

// #include <deal.II/fe/fe_q.h>
// #include <deal.II/fe/fe_system.h>
// #include <deal.II/fe/fe_values.h>
// #include <deal.II/fe/mapping_q1.h>
// #include <deal.II/fe/mapping_q1_eulerian.h>

// #include <deal.II/numerics/data_out.h>
// #include <deal.II/numerics/solution_transfer.h>
// #include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "bem_fma.h"
#include "bem_problem.h"
#include "boundary_conditions.h"
#include "computational_domain.h"

#include <sys/time.h>
#include "Teuchos_TimeMonitor.hpp"

using namespace dealii;

template <int dim>
class Driver : public ParameterAcceptor
{
public:
  Driver();

  ~Driver();

  /// method to declare the parameters
  /// to be read from the parameters file

  virtual void
  declare_parameters(ParameterHandler &prm);

  /// method to parse the needed parameters
  /// from the parameters file

  virtual void
  parse_parameters(ParameterHandler &prm);


  void
  run();

protected:
  ConditionalOStream pcout;

  MPI_Comm mpi_communicator;

  ComputationalDomain<dim> computational_domain;

  BEMProblem<dim> bem_problem;

  BoundaryConditions<dim> boundary_conditions;

  ParameterHandler prm;

  bool global_refinement;

  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;
};

using Teuchos::RCP;
using Teuchos::Time;
using Teuchos::TimeMonitor;

RCP<Time> TotalTime  = Teuchos::TimeMonitor::getNewTimer("Total Time");
RCP<Time> MeshTime   = Teuchos::TimeMonitor::getNewTimer("Mesh Time");
RCP<Time> OutputTime = Teuchos::TimeMonitor::getNewTimer("Output Time");
RCP<Time> SolveTime  = Teuchos::TimeMonitor::getNewTimer("Solve Time");

using namespace std;

template <int dim>
Driver<dim>::Driver()
  : pcout(std::cout)
  , mpi_communicator(MPI_COMM_WORLD)
  , computational_domain(mpi_communicator)
  , bem_problem(computational_domain, mpi_communicator)
  , boundary_conditions(computational_domain, bem_problem)
  , prm()
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
{
  pcout.set_condition(this_mpi_process == 0);
}

template <int dim>
Driver<dim>::~Driver()
{}


template <int dim>
void
Driver<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.declare_entry("Set Global Refinement", "true", Patterns::Bool());
}

template <int dim>
void
Driver<dim>::parse_parameters(ParameterHandler &prm)
{
  global_refinement = prm.get_bool("Set Global Refinement");
}


template <int dim>
void
Driver<dim>::run()
{
  Teuchos::TimeMonitor LocalTimer(*TotalTime);

  computational_domain.read_domain();
  computational_domain.refine_and_resize(computational_domain.n_cycles);
  computational_domain.update_triangulation();

  bem_problem.reinit();
  boundary_conditions.initialize();

  boundary_conditions.assign_potential();
  boundary_conditions.assign_potential_normal_derivative();

  bem_problem.compute_hypersingular_free_coeffs();
  boundary_conditions.compute_errors();
  boundary_conditions.output_results(boundary_conditions.output_file_name);

  Teuchos::TimeMonitor::summarize();
}

template class Driver<2>;
template class Driver<3>;


int main(int argc, char *argv[])
{
  try
    {
      unsigned int threads;
      if (argc == 1)
        threads = numbers::invalid_unsigned_int;
      else
        threads = atoi(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, threads);

      std::string pname =
        "parameters_bem_" + std::to_string(DEAL_II_DIMENSION) + ".prm";
      std::string pname2 =
        "used_parameters_bem_" + std::to_string(DEAL_II_DIMENSION) + ".prm";

      Driver<DEAL_II_DIMENSION> driver;
      dealii::ParameterAcceptor::initialize(pname, pname2);

      driver.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
