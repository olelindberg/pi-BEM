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

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "./../tests.h"
#include "bem_problem.h"
#include "boundary_conditions.h"
#include "computational_domain.h"
#include "surface_integral_util.h"

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
  Options options(argc, argv);
  std::cout << "Input file: " << options.get_input_file() << "\n";

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  // initlog();
  MPI_Comm                 mpi_communicator(MPI_COMM_WORLD);
  ComputationalDomain<3>   computational_domain(mpi_communicator);
  surface_integral_util<3> surf_integr(computational_domain, mpi_communicator);

  deal2lkit::ParameterAcceptor::initialize(options.get_input_file(), "used.prm");

  computational_domain.read_domain();
  computational_domain.read_cad_files();
  computational_domain.assign_manifold_projectors(1.0e-6);
  computational_domain.getTria().refine_global(2);

  surf_integr.reinit();
  std::vector<double> areas = surf_integr.ssurffint();
  for (auto &area : areas)
    std::cout << "area,               A  = " << area << std::endl;
}
