#include "../include/driver.h"

#include <boost/filesystem.hpp>

#include <sys/time.h>

#include <fstream>

#include "../include/surface_integral_util.h"
#include "JSON_BodyReader.h"
#include "Teuchos_TimeMonitor.hpp"

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
Driver<dim>::run(std::string input_path, std::string output_path)
{
  // Constant parameters:
  const double gravity = 9.80665;
  const double density = 1000;

  // Read and create body/ship:
  std::string bodyfilename = boost::filesystem::path(input_path).append("body.json").string();
  Body        body;
  if (!JSON_BodyReader::read(bodyfilename, body))
  {
    std::cout << "Reading json body file failed ..." << std::endl;
    return;
  }
  body.print();


  {
    Teuchos::TimeMonitor LocalTimer(*TotalTime);
    unsigned int         local_refinement_cycles = 0;
    {
      computational_domain.read_domain(input_path);
      if (global_refinement)
      {
        std::cout << "Global refinement ...\n";
        computational_domain.refine_and_resize(computational_domain.n_cycles, input_path);
      }
      else
      {
        std::cout << "---------------------------------------------------------\n";
        std::cout << "Pre adaptive global refinement ...\n";
        computational_domain.refine_and_resize(computational_domain.pre_global_refinements,
                                               input_path);
        local_refinement_cycles = computational_domain.n_cycles;
      }
    }

    computational_domain.update_triangulation();
    for (unsigned int i = 0; i <= local_refinement_cycles; ++i)
    {
      std::cout << "-----------------------------------------------------------\n";
      std::cout << "Adaptive refinement " << i << " ...\n";
      {
        Teuchos::TimeMonitor LocalTimer(*SolveTime);
        bem_problem.reinit();
        boundary_conditions.solve_problem(body);
      }
      if (!global_refinement && i < local_refinement_cycles)
      {
        // Compute error estimator and local refinement strategy
        bem_problem.adaptive_refinement(boundary_conditions.get_phi());
        computational_domain.update_triangulation();
      }
    }


    //-------------------------------------------------------------------------
    // Waterplane area and displacement:
    //-------------------------------------------------------------------------
    auto area   = bem_problem.area_integral(body);
    auto volume = bem_problem.volume_integral(body);

    //-------------------------------------------------------------------------
    // Hydrostatic and hydrodynamic pressures:
    //-------------------------------------------------------------------------
    bem_problem.hydrostatic_pressure(gravity,
                                     density,
                                     body.getDraft(),
                                     boundary_conditions.get_hydrostatic_pressure());

    bem_problem.hydrodynamic_pressure(density,
                                      boundary_conditions.get_wind(),
                                      boundary_conditions.get_hydrodynamic_pressure());


    //-------------------------------------------------------------------------
    // Pressure centers:
    //-------------------------------------------------------------------------
    Tensor<1, dim, double> hydrostatic_pressure_center;
    bem_problem.center_of_pressure(body,
                                   boundary_conditions.get_hydrostatic_pressure(),
                                   hydrostatic_pressure_center);

    Tensor<1, dim, double> hydrodynamic_pressure_center;
    bem_problem.center_of_pressure(body,
                                   boundary_conditions.get_hydrodynamic_pressure(),
                                   hydrodynamic_pressure_center);

    std::cout << "Waterplane area        : A   = " << area << "\n";
    std::cout << "Displacement           : V   = " << volume << "\n";
    std::cout << "Center of flotation    : COF = " << hydrostatic_pressure_center << "\n";
    std::cout << "Dynamic pressure center: COP = " << hydrodynamic_pressure_center << "\n";

    //-------------------------------------------------------------------------
    // Pressure Forces:
    //-------------------------------------------------------------------------
    Tensor<1, dim, double> hydrostaticForce;
    Tensor<1, dim, double> hydrostaticMoment;
    bem_problem.pressure_force_and_moment(body,
                                          boundary_conditions.get_hydrostatic_pressure(),
                                          hydrostaticForce,
                                          hydrostaticMoment);

    Tensor<1, dim, double> hydrodynamicForce;
    Tensor<1, dim, double> hydrodynamicMoment;
    bem_problem.pressure_force_and_moment(body,
                                          boundary_conditions.get_hydrodynamic_pressure(),
                                          hydrodynamicForce,
                                          hydrodynamicMoment);


    std::vector<Point<dim>> elevation;
    bem_problem.free_surface_elevation(
      gravity, density, body, boundary_conditions.get_hydrodynamic_pressure(), elevation);

    //-------------------------------------------------------------------------
    // Save forces:
    //-------------------------------------------------------------------------
    std::fstream file;
    std::string  force_filename = boost::filesystem::path(output_path).append("force.csv").string();
    file.open(force_filename, std::fstream::out);
    if (file.is_open())
    {
      file << "# Fx [N], Fy [N], Fz [N], Mx [Nm], My [Nm], Mz [Nm]\n";
      file << hydrodynamicForce[0] << ", " << hydrodynamicForce[1] << ", " << hydrodynamicForce[2]
           << ", " << hydrodynamicMoment[0] << ", " << hydrodynamicMoment[1] << ", "
           << hydrodynamicMoment[2] << "\n";
      file << hydrostaticForce[0] << ", " << hydrostaticForce[1] << ", " << hydrostaticForce[2]
           << ", " << hydrostaticMoment[0] << ", " << hydrostaticMoment[1] << ", "
           << hydrostaticMoment[2] << "\n";
      file.close();
    }

    //-------------------------------------------------------------------------
    // Save wave elevation:
    //-------------------------------------------------------------------------
    {
      std::fstream file;
      std::string  filename = boost::filesystem::path(output_path).append("elevation.csv").string();
      file.open(filename, std::fstream::out);
      if (file.is_open())
      {
        file << "# x [m], y [m], z [m]\n";
        for (auto &elev : elevation)
          file << elev[0] << ", " << elev[1] << ", " << elev[2] << "\n";
        file.close();
      }
    }

    std::string filename =
      boost::filesystem::path(output_path).append(boundary_conditions.output_file_name).string();
    boundary_conditions.output_results(filename);
  }

  // Write a summary of all timers
  //  Teuchos::TimeMonitor::summarize();
}

// template class Driver<2>;
template class Driver<3>;
