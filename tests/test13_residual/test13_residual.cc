//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Read a file in iges format, and write it out again in the same
// format.

#include "./../tests.h"
#include "bem_problem.h"
#include "boundary_conditions.h"
#include "computational_domain.h"


#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>
#include <boost/program_options.hpp>
#include <filesystem>


class Options
{
public:
  Options(int argc, char *argv[])
  {
    _desc = std::make_shared<boost::program_options::options_description>(
      "test08_center_of_pressure options");
    _desc->add_options()("help", "produce help message");
    _desc->add_options()("input-file", boost::program_options::value<std::string>(), "Input file");


    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, *_desc),
                                  _vm);
    boost::program_options::notify(_vm);

    if (_vm.count("help"))
      std::cout << *_desc << "\n";
  }

  std::string
  get_input_file()
  {
    return _vm["input-file"].as<std::string>();
  };

private:
  std::shared_ptr<boost::program_options::options_description> _desc;
  boost::program_options::variables_map                        _vm;
};


int
main(int argc, char **argv)
{
  std::string thisFileName {__FILE__};
  std::filesystem::path cwd = std::filesystem::path(thisFileName).parent_path();
  
 
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  // initlog();
  MPI_Comm               mpi_communicator(MPI_COMM_WORLD);
  auto  computational_domain = ComputationalDomain<3>::create(mpi_communicator);
  BEMProblem<3>          bem_problem(computational_domain, mpi_communicator);
  BoundaryConditions<3>  boundary_conditions(computational_domain, bem_problem);

  deal2lkit::ParameterAcceptor::initialize(cwd.append("../../NACA_FILES/parameters_bem_3.prm"), "used.prm");

  computational_domain->read_domain();
  computational_domain->read_cad_files();
  computational_domain->assign_manifold_projectors(1.0e-6);
  computational_domain->getTria().refine_global(2);

  bem_problem.reinit();

 // bem_problem.dynamic_pressure(boundary_conditions.get_wind(), boundary_conditions.get_pressure());

//  std::vector<double>               areas   = bem_problem.area_integral();
////  std::vector<Tensor<1, 3, double>> volumes = bem_problem.volume_integral();

//  std::vector<Tensor<1, 3, double>> pressure_centers;
//  bem_problem.center_of_pressure(boundary_conditions.get_pressure(), pressure_centers);

//  for (auto &area : areas)
//    std::cout << "area,                 A  = " << area << std::endl;
//  for (auto &volume : volumes)
//    std::cout << "volume,               V  = " << volume << std::endl;
//  for (auto &cp : pressure_centers)
//    std::cout << "center of pressure,   cp = " << cp << std::endl;

  //  std::string   filename0 = ("meshResult.inp");
  //  std::ofstream logfile0(filename0.c_str());
  //  GridOut       grid_out0;
  //  grid_out0.write_ucd(computational_domain.getTria(), logfile0);

  
}
