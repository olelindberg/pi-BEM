#include "../include/driver.h"

#include <boost/filesystem.hpp>

#include <sys/time.h>

#include <fstream>

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
        boundary_conditions.solve_problem();
      }
      if (!global_refinement && i < local_refinement_cycles)
      {
        // Compute error estimator and local refinement strategy
        bem_problem.adaptive_refinement(boundary_conditions.get_phi());
        computational_domain.update_triangulation();
      }
    }

    bem_problem.dynamic_pressure(boundary_conditions.get_wind(),
                                 boundary_conditions.get_pressure());

    std::vector<Point<dim>> points;
    std::vector<Point<dim>> velocities;

    double Lx = 4.0;
    double Ly = 4.0;
    int    Nx = 21;
    int    Ny = 21;

    double dx = 1.0;
    if (Nx > 1)
      dx = Lx / (Nx - 1);

    double dy = 1.0;
    if (Ny > 1)
      dy = Ly / (Ny - 1);

    double x0 = -2.0;
    double y0 = -2.0;
    double z0 = 0.0;
    for (int j = 0; j < Ny; ++j)
    {
      for (int i = 0; i < Nx; ++i)
      {
        double     x = x0 + dx * i;
        double     y = y0 + dy * j;
        double     z = z0;
        Point<dim> pnt(x, y, z);
        points.push_back(pnt);
      }
    }

    std::vector<double> potentials;
    bem_problem.velocity(0.0,
                         boundary_conditions.get_phi(),
                         boundary_conditions.get_dphi_dn(),
                         points,
                         potentials,
                         velocities);

    auto                    area   = bem_problem.area_integral();
    auto                    volume = bem_problem.volume_integral();
    auto                    forces = bem_problem.pressure_force(boundary_conditions.get_pressure());
    std::vector<Point<dim>> elevation;
    bem_problem.free_surface_elevation(boundary_conditions.get_pressure(), elevation);

    std::cout << "area  : " << area << "\n";
    std::cout << "volume: " << volume[0] << ", " << volume[1] << ", " << volume[2] << "\n";

    //-------------------------------------------------------------------------
    // Save forces:
    //-------------------------------------------------------------------------
    {
      std::fstream file;
      std::string  filename = boost::filesystem::path(output_path).append("velocity.csv").string();
      file.open(filename, std::fstream::out);
      if (file.is_open())
      {
        file << "# x [m], y [m], z [m], phi [m^2/s], u [m/s], v [m/s], w [m/s] \n";
        for (int i = 0; i < points.size(); ++i)
        {
          file << points[i][0] << ", " << points[i][1] << ", " << points[i][2] << ", "
               << potentials[i] << ", " << velocities[i][0] << ", " << velocities[i][1] << ", "
               << velocities[i][2] << "\n";
        }
        file.close();
      }
    }
    //-------------------------------------------------------------------------
    // Save forces:
    //-------------------------------------------------------------------------
    std::fstream file;
    std::string force_filename = boost::filesystem::path(output_path).append("forces.csv").string();
    file.open(force_filename, std::fstream::out);
    if (file.is_open())
    {
      file << "# material id, Fx [N], Fy [N], Fz [N]\n";
      int cnt = 0;
      for (auto &force : forces)
      {
        file << cnt << ", " << force[0] << ", " << force[1] << ", " << force[2] << "\n";
        ++cnt;
      }
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

    boundary_conditions.compute_errors(output_path);

    std::string filename =
      boost::filesystem::path(output_path).append(boundary_conditions.output_file_name).string();
    boundary_conditions.output_results(filename);
  }

  // Write a summary of all timers
  //  Teuchos::TimeMonitor::summarize();
}

template class Driver<2>;
template class Driver<3>;
