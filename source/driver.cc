#include "../include/driver.h"

#include <boost/filesystem.hpp>

#include <sys/time.h>

#include <fstream>

#include "../include/surface_integral_util.h"
#include "JSON_BodyReader.h"
#include "Teuchos_TimeMonitor.hpp"

#include "../include/WireUtil.h"
#include <BRepBuilderAPI_MakeWire.hxx>

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
  //-------------------------------------------------------------------------
  // Pre steps:
  //-------------------------------------------------------------------------

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

  if (global_refinement)
  {
    std::cout << "Global refinement ...\n";
    Teuchos::TimeMonitor LocalTimer(*TotalTime);
    computational_domain.read_domain(input_path);
    computational_domain.refine_and_resize(input_path);
    computational_domain.update_triangulation();
    bem_problem.reinit();
    boundary_conditions.solve_problem(body);
  }
  else // adaptive refinement
  {
    pcout << "Adaptive refinement ...\n";
    Teuchos::TimeMonitor LocalTimer(*TotalTime);
    computational_domain.read_domain(input_path);
    computational_domain.refine_and_resize(input_path);
    computational_domain.update_triangulation();
    bem_problem.reinit();
    boundary_conditions.solve_problem(body);

    for (unsigned int i = 0; i < computational_domain.n_cycles; ++i)
    {
      pcout << "Refinement level " << i << " ...\n";

      bem_problem.hydrodynamic_pressure(density,
                                        boundary_conditions.get_wind(),
                                        boundary_conditions.get_hydrodynamic_pressure());

      bem_problem.adaptive_refinement(boundary_conditions.get_hydrodynamic_pressure());
      computational_domain.update_triangulation();

      bem_problem.reinit();
      boundary_conditions.solve_problem(body);
    }
  }


  //-------------------------------------------------------------------------
  // Main steps:
  //-------------------------------------------------------------------------
  {
    //-------------------------------------------------------------------------
    // Post steps:
    //-------------------------------------------------------------------------


    //-------------------------------------------------------------------------
    // Waterplane area and displacement:
    //-------------------------------------------------------------------------
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

    std::cout << "Displacement           : V = " << volume << "\n";
    std::cout << "Static pressure center :     " << hydrostatic_pressure_center << "\n";
    std::cout << "Dynamic pressure center:     " << hydrodynamic_pressure_center << "\n";

    if (false)
    {
      BRepBuilderAPI_MakeWire wirebuilder;
      for (auto id : body.getWaterlineIndices())
        wirebuilder.Add(TopoDS::Wire(computational_domain.cad_curves[id - 11]));
      if (wirebuilder.IsDone())
      {
        double         x0 = body.getCenterOfGravity()[0];
        double         y0 = body.getCenterOfGravity()[1];
        SurfaceMoments sm(x0, y0);
        if (!WireUtil::surfaceMoments(wirebuilder.Wire(), sm))
          pcout << "Surface moments failed ..." << std::endl;

        pcout << "x0  : " << sm.getx0() << std::endl;
        pcout << "y0  : " << sm.gety0() << std::endl;
        pcout << "S0  : " << sm.getS0() << std::endl;
        pcout << "Sx  : " << sm.getSx() << std::endl;
        pcout << "Sy  : " << sm.getSy() << std::endl;
        pcout << "Sxx : " << sm.getSxx() << std::endl;
        pcout << "Sxy : " << sm.getSxy() << std::endl;
        pcout << "Syy : " << sm.getSyy() << std::endl;
      }
    }
    // dealii::Point<3> tmp;
    // for (int i = 0; i < 3; ++i)
    //   tmp[i] = hydrostatic_pressure_center[i];
    // WaterPlaneMoments wpm = bem_problem.water_plane_moments(body, tmp);

    // std::cout << "S0  = " << wpm.getS0() << "\n";
    // std::cout << "Sx  = " << wpm.getSx() << "\n";
    // std::cout << "Sy  = " << wpm.getSy() << "\n";
    // std::cout << "Sxx = " << wpm.getSxx() << "\n";
    // std::cout << "Sxy = " << wpm.getSxy() << "\n";
    // std::cout << "Syy = " << wpm.getSyy() << "\n";



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
