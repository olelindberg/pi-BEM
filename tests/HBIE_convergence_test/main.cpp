#include <sys/time.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "Teuchos_TimeMonitor.hpp"
#include "bem_fma.h"
#include "bem_problem.h"
#include "boundary_conditions.h"
#include "computational_domain.h"

using namespace dealii;

class DriverSetup
{
public:
  DriverSetup(int val1, int val2, int val3, int val4, int val5)
    : n_cycles(val1)
    , fe_order(val2)
    , quadrature_order(val3)
    , rho_quadrature_order(val4)
    , theta_quadrature_order(val5)
  {}
  unsigned int n_cycles               = 2;
  unsigned int fe_order               = 2;
  unsigned int quadrature_order       = 2;
  unsigned int rho_quadrature_order   = 2;
  unsigned int theta_quadrature_order = 2;
};

template <int dim>
class Driver : public ParameterAcceptor
{
public:
  Driver(DriverSetup setup);

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

  ConditionalOStream pcout;

  MPI_Comm mpi_communicator;

  ComputationalDomain<dim> computational_domain;

  BEMProblem<dim> bem_problem;

  BoundaryConditions<dim> boundary_conditions;

  ParameterHandler prm;

  bool global_refinement;

  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  double
  get_hbietime()
  {
    return _hbietime;
  }

private:
  DriverSetup _setup;
  double      _hbietime;
};

using Teuchos::RCP;
using Teuchos::Time;
using Teuchos::TimeMonitor;

RCP<Time> TotalTime  = Teuchos::TimeMonitor::getNewTimer("Total Time");
RCP<Time> HBIETime   = Teuchos::TimeMonitor::getNewTimer("HBIE Time");
RCP<Time> MeshTime   = Teuchos::TimeMonitor::getNewTimer("Mesh Time");
RCP<Time> OutputTime = Teuchos::TimeMonitor::getNewTimer("Output Time");
RCP<Time> SolveTime  = Teuchos::TimeMonitor::getNewTimer("Solve Time");

using namespace std;

template <int dim>
Driver<dim>::Driver(DriverSetup setup)
  : pcout(std::cout)
  , mpi_communicator(MPI_COMM_WORLD)
  , computational_domain(mpi_communicator)
  , bem_problem(computational_domain, mpi_communicator)
  , boundary_conditions(computational_domain, bem_problem)
  , prm()
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , _setup(setup)
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
  computational_domain.refine_and_resize(_setup.n_cycles);
  computational_domain.update_triangulation();

  bem_problem.set_mapping_degree(_setup.fe_order);
  bem_problem.set_scalar_fe_order(_setup.fe_order);
  bem_problem.set_vector_fe_order(_setup.fe_order);
  bem_problem.set_quadrature_order(_setup.quadrature_order);
  bem_problem.set_hbie_radial_quadrature_order(_setup.rho_quadrature_order);
  bem_problem.set_hbie_angular_quadrature_order(_setup.theta_quadrature_order);

  bem_problem.reinit();
  boundary_conditions.initialize();

  boundary_conditions.assign_potential();
  boundary_conditions.assign_potential_normal_derivative();

  HBIETime->enable();
  HBIETime->start(true);
  bem_problem.compute_hypersingular_free_coeffs();
  bem_problem.compute_gradients_hypersingular(boundary_conditions.get_phi(), boundary_conditions.get_dphi_dn());
  _hbietime = HBIETime->stop();

  boundary_conditions.compute_errors();
  boundary_conditions.output_results(boundary_conditions.output_file_name);

  Teuchos::TimeMonitor::summarize();
}

template class Driver<2>;
template class Driver<3>;


int
main(int argc, char *argv[])
{
  try
    {
      unsigned int threads;
      if (argc == 1)
        threads = numbers::invalid_unsigned_int;
      else
        threads = atoi(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, threads);

      std::string pname  = "parameters_bem_" + std::to_string(DEAL_II_DIMENSION) + ".prm";
      std::string pname2 = "used_parameters_bem_" + std::to_string(DEAL_II_DIMENSION) + ".prm";

      std::vector<DriverSetup> setups;
      setups.push_back(DriverSetup(2, 2, 5, 4, 20));
      // setups.push_back(DriverSetup(0, 1, 20, 4, 20));
      // setups.push_back(DriverSetup(1, 1, 20, 4, 20));
      // setups.push_back(DriverSetup(2, 1, 20, 4, 20));
      // setups.push_back(DriverSetup(3, 1, 20, 4, 20));
      // setups.push_back(DriverSetup(4, 1, 20, 4, 20));
      // setups.push_back(DriverSetup(0, 2, 20, 4, 20));
      // setups.push_back(DriverSetup(1, 2, 20, 4, 20));
      // setups.push_back(DriverSetup(2, 2, 20, 4, 20));
      // setups.push_back(DriverSetup(3, 2, 20, 4, 20));
      // setups.push_back(DriverSetup(4, 2, 20, 4, 20));
      // setups.push_back(DriverSetup(0, 3, 20, 4, 20));
      // setups.push_back(DriverSetup(1, 3, 20, 4, 20));
      // setups.push_back(DriverSetup(2, 3, 20, 4, 20));
      // setups.push_back(DriverSetup(3, 3, 20, 4, 20));
      // setups.push_back(DriverSetup(0, 4, 20, 4, 20));
      // setups.push_back(DriverSetup(1, 4, 20, 4, 20));
      // setups.push_back(DriverSetup(2, 4, 20, 4, 20));


      std::ofstream file;
      file.open("HBIE_convergence_test.csv", std::ios::trunc);
      file << "# ncells, ndofs, degree, quadrature_order, rho_quadrature_order, theta_quadrature_order, grad_L2_error, grad_inf_error, hbie_time \n";
      file.close();
      for (auto setup : setups)
        {
          Driver<DEAL_II_DIMENSION> driver(setup);
          dealii::ParameterAcceptor::initialize(pname, pname2);

          driver.run();

          std::ofstream file;
          file.open("HBIE_convergence_test.csv", std::ios::app);
          file << driver.computational_domain.tria.n_active_cells() << ", " << driver.bem_problem.dh.n_dofs() << ", " << driver.bem_problem.fe->degree << ", " << setup.quadrature_order << ", " << setup.rho_quadrature_order << ", " << setup.theta_quadrature_order << ", " << driver.boundary_conditions.get_grad_phi_L2_error() << ", " << driver.boundary_conditions.get_grad_phi_max_error() << ", " << driver.get_hbietime() << "\n";
          file.close();
        }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl << "----------------------------------------------------" << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what() << std::endl << "Aborting!" << std::endl << "----------------------------------------------------" << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl << "----------------------------------------------------" << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!" << std::endl << "----------------------------------------------------" << std::endl;
      return 1;
    }

  return 0;
}
