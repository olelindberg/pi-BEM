#include "../include/driver.h"
#include "Teuchos_TimeMonitor.hpp"
#include <boost/filesystem.hpp>
#include <sys/time.h>

#include <fstream>

using Teuchos::RCP;
using Teuchos::Time;
using Teuchos::TimeMonitor;

RCP<Time> TotalTime  = Teuchos::TimeMonitor::getNewTimer ("Total Time");
RCP<Time> MeshTime   = Teuchos::TimeMonitor::getNewTimer ("Mesh Time");
RCP<Time> OutputTime = Teuchos::TimeMonitor::getNewTimer ("Output Time");
RCP<Time> SolveTime  = Teuchos::TimeMonitor::getNewTimer ("Solve Time");

using namespace std;

template <int dim>
Driver<dim>::Driver ()
    : pcout (std::cout), mpi_communicator (MPI_COMM_WORLD), computational_domain (mpi_communicator), bem_problem (computational_domain, mpi_communicator),
      boundary_conditions (computational_domain, bem_problem), prm (), n_mpi_processes (Utilities::MPI::n_mpi_processes (mpi_communicator)),
      this_mpi_process (Utilities::MPI::this_mpi_process (mpi_communicator))
{
  pcout.set_condition (this_mpi_process == 0);
}

template <int dim> Driver<dim>::~Driver () {}

template <int dim>
void
Driver<dim>::declare_parameters (ParameterHandler &prm)
{
  prm.declare_entry ("Set Global Refinement", "true", Patterns::Bool ());
}

template <int dim>
void
Driver<dim>::parse_parameters (ParameterHandler &prm)
{
  global_refinement = prm.get_bool ("Set Global Refinement");
}

template <int dim>
void
Driver<dim>::run (std::string input_path, std::string output_path)
{

    std::vector<int> force_material_ids = {1,2,3,4};

  {
    Teuchos::TimeMonitor LocalTimer (*TotalTime);
    unsigned int         local_refinement_cycles = 0;
    {
      computational_domain.read_domain (input_path);
      if (global_refinement)
      {
        std::cout << "Global refinement ...\n";
        computational_domain.refine_and_resize (computational_domain.n_cycles, input_path);
      }
      else
      {
        std::cout << "---------------------------------------------------------\n";
        std::cout << "Pre adaptive global refinement ...\n";
        computational_domain.refine_and_resize (computational_domain.pre_global_refinements, input_path);
        local_refinement_cycles = computational_domain.n_cycles;
      }
    }

    computational_domain.update_triangulation ();
    for (unsigned int i = 0; i <= local_refinement_cycles; ++i)
    {
      std::cout << "-----------------------------------------------------------\n";
      std::cout << "Adaptive refinement " << i << " ...\n";
      {
        Teuchos::TimeMonitor LocalTimer (*SolveTime);
        bem_problem.reinit ();
        boundary_conditions.solve_problem ();
      }
      if (!global_refinement && i < local_refinement_cycles)
      {
        // Compute error estimator and local refinement strategy
        bem_problem.adaptive_refinement (boundary_conditions.get_phi ());
        computational_domain.update_triangulation ();
      }
    }

    bem_problem.dynamic_pressure (boundary_conditions.get_wind (), boundary_conditions.get_pressure ());


    auto force = bem_problem.pressure_force (boundary_conditions.get_pressure (), force_material_ids);

    
//    std::to_string(force_material_id);
    std::fstream file;
    std::string force_filename = boost::filesystem::path (output_path).append ("force.csv").string ();
    file.open (force_filename, std::fstream::out);
    if (file.is_open ())
    {
      file << "# Fx [N], Fy [N], Fz [N]\n";
      file << force[0] << ", " << force[1] << ", " << force[2];
      file.close ();
    }

    boundary_conditions.compute_errors (output_path);

    std::string filename = boost::filesystem::path (output_path).append (boundary_conditions.output_file_name).string ();
    boundary_conditions.output_results (filename);
  }

  // Write a summary of all timers
  Teuchos::TimeMonitor::summarize ();
}

template class Driver<2>;
template class Driver<3>;
