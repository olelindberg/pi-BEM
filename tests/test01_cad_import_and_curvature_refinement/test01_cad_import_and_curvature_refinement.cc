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

#include "computational_domain.h"
#include "./../tests.h"

#include <filesystem>

int
main(int argc, char **argv)
{
  std::string thisFileName {__FILE__};
  std::filesystem::path cwd = std::filesystem::path(thisFileName).parent_path();
  std::filesystem::current_path(cwd);

  std::filesystem::path dataPath = std::filesystem::path("../test_data/NACA_FILES/");


  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  // initlog();
  MPI_Comm               mpi_communicator(MPI_COMM_WORLD);
  ComputationalDomain<3> computational_domain(mpi_communicator);

  deal2lkit::ParameterAcceptor::initialize(cwd / std::filesystem::path("parameters_bem_3.prm"), "used.prm");  
  computational_domain.read_domain(dataPath);
  computational_domain.read_cad_files(dataPath);
  computational_domain.refine_and_resize(cwd);
}
