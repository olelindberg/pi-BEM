#include "../include/driver.h"

#include <boost/filesystem.hpp>

#include <sys/time.h>

#include <fstream>

#include "../include/AdaptiveRefinement.h"
#include "../include/Writer.h"
#include "../include/surface_integral_util.h"

#include "BodySettingsReaderJSON.h"
#include "PiBEMSettingsReaderJSON.h"

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
  try
  {
    //-------------------------------------------------------------------------
    // Pre steps:
    //-------------------------------------------------------------------------
    std::string   pibemFilename = boost::filesystem::path(input_path).append("pibem.json").string();
    PiBEMSettings pibemSettings;
    if (!PiBEMSettingsReaderJSON::read(pibemFilename, pibemSettings))
    {
      std::cout << "Reading json pibem settings file failed ..." << std::endl;
      return;
    }
    pibemSettings.print();

    // Read and create body/ship:
    std::string bodyfilename = boost::filesystem::path(input_path).append("body.json").string();
    Body        body;
    if (!BodySettingsReaderJSON::read(bodyfilename, body))
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

      AdaptiveRefinement adaptiveRefinement(pcout,
                                            mpi_communicator,
                                            pibemSettings.potentialErrorEstimatorMax,
                                            pibemSettings.velocityErrorEstimatorMax,
                                            pibemSettings.aspectRatioMax,
                                            pibemSettings.cellSizeMin,
                                            pibemSettings.iterMax);

      for (int i = 0; i < pibemSettings.adaptiveRefinementLevels; ++i)
      {
        pcout << "Refinement level " << i << " ...\n";



        if (!adaptiveRefinement.refine(n_mpi_processes,
                                       bem_problem.this_mpi_process,
                                       *bem_problem.fe,
                                       *bem_problem.gradient_fe,
                                       bem_problem.dh,
                                       bem_problem.gradient_dh,
                                       boundary_conditions.get_phi(),
                                       bem_problem.vector_gradients_solution,
                                       computational_domain.tria))
        {
          break;
        }
        pcout << "Degrees of Freedom: DOF = " << bem_problem.dh.n_dofs() << std::endl;



        if (this_mpi_process == 0)
        {
          std::fstream file;
          std::string  filename = boost::filesystem::path(output_path)
                                   .append(std::string("potentialErrorEstimatorLevel")
                                             .append(std::to_string(i))
                                             .append(".csv"))
                                   .string();
          file.open(filename, std::fstream::out);
          if (file.is_open())
          {
            for (auto &val : adaptiveRefinement.get_error_estimator_potential())
              file << val << "\n";
            file.close();
          }
        }
        pcout << "Degrees of Freedom: DOF = " << bem_problem.dh.n_dofs() << std::endl;

        if (this_mpi_process == 0)
        {
          std::fstream file;
          std::string  filename =
            boost::filesystem::path(output_path)
              .append(
                std::string("velocityErrorEstimatorLevel").append(std::to_string(i)).append(".csv"))
              .string();
          file.open(filename, std::fstream::out);
          if (file.is_open())
          {
            for (auto &val : adaptiveRefinement.get_error_estimator_velocity())
              file << val << "\n";
            file.close();
          }
        }
        pcout << "Degrees of Freedom: DOF = " << bem_problem.dh.n_dofs() << std::endl;

        computational_domain.update_triangulation();
        pcout << "Degrees of Freedom: DOF = " << bem_problem.dh.n_dofs() << std::endl;

        bem_problem.reinit();
        pcout << "Degrees of Freedom: DOF = " << bem_problem.dh.n_dofs() << std::endl;

        Writer writer;
        writer.addScalarField("error_estimator",
                              adaptiveRefinement.get_error_estimator_potential());
        writer.saveScalarFields(std::string(input_path).append("/scalars.vtu"),
                                bem_problem.dh,
                                bem_problem.mapping,
                                bem_problem.mapping_degree);

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
      bem_problem.hydrostatic_pressure(pibemSettings.gravity,
                                       pibemSettings.density,
                                       body.getDraft(),
                                       boundary_conditions.get_hydrostatic_pressure());

      bem_problem.hydrodynamic_pressure(pibemSettings.density,
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
      bem_problem.free_surface_elevation(pibemSettings.gravity,
                                         pibemSettings.density,
                                         body,
                                         boundary_conditions.get_hydrodynamic_pressure(),
                                         elevation);

      if (this_mpi_process == 0)
      {
        //-------------------------------------------------------------------------
        // Save forces:
        //-------------------------------------------------------------------------
        std::fstream file;
        std::string  force_filename =
          boost::filesystem::path(output_path).append("force.csv").string();
        file.open(force_filename, std::fstream::out);
        if (file.is_open())
        {
          file << "# Fx [N], Fy [N], Fz [N], Mx [Nm], My [Nm], Mz [Nm]\n";
          file << hydrodynamicForce[0] << ", " << hydrodynamicForce[1] << ", "
               << hydrodynamicForce[2] << ", " << hydrodynamicMoment[0] << ", "
               << hydrodynamicMoment[1] << ", " << hydrodynamicMoment[2] << "\n";
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
          std::string  filename =
            boost::filesystem::path(output_path).append("elevation.csv").string();
          file.open(filename, std::fstream::out);
          if (file.is_open())
          {
            file << "# x [m], y [m], z [m]\n";
            for (auto &elev : elevation)
              file << elev[0] << ", " << elev[1] << ", " << elev[2] << "\n";
            file.close();
          }
        }
      }
      std::string filename =
        boost::filesystem::path(output_path).append(boundary_conditions.output_file_name).string();
      boundary_conditions.output_results(filename); // \todo change to Writer in Writer.h
    }
  }
  catch (std::exception &e)
  {
    std::string filename =
      boost::filesystem::path(output_path)
        .append(std::string("error_dump_").append(boundary_conditions.output_file_name))
        .string();
    boundary_conditions.output_results(filename); // \todo change to Writer in Writer.h
    std::cout << e.what() << std::endl;
    return;
  }
  // Write a summary of all timers
  //  Teuchos::TimeMonitor::summarize();
}

// template class Driver<2>;
template class Driver<3>;
