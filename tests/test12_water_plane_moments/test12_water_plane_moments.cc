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

#include "./../tests.h"
#include "bem_problem.h"
#include "computational_domain.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  // initlog();
  MPI_Comm               mpi_communicator(MPI_COMM_WORLD);
  ComputationalDomain<3> computational_domain(mpi_communicator);
  BEMProblem<3>          bem_problem(computational_domain, mpi_communicator);

  deal2lkit::ParameterAcceptor::initialize("test06_area_and_volume.prm", "used.prm");

  computational_domain.read_domain();
  computational_domain.getTria().refine_global(3);

  bem_problem.reinit();

  Body             body;
  std::vector<int> wl = {11, 12};
  body.setWaterlineIndices(wl);

  auto area   = bem_problem.area_integral(body);
  auto volume = bem_problem.volume_integral(body);

  std::cout << "area   = " << area << std::endl;
  std::cout << "volume = " << volume << std::endl;

  std::string   filename0 = ("meshResult.inp");
  std::ofstream logfile0(filename0.c_str());
  GridOut       grid_out0;
  grid_out0.write_ucd(computational_domain.getTria(), logfile0);
}
