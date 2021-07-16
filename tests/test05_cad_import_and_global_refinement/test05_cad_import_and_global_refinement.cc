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
#include <deal.II/grid/grid_tools.h>

enum class REFINE {GLOBAL,ASPECT_RATIO};

int
main(int argc, char **argv)
{
  auto refine = REFINE::GLOBAL;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  // initlog();
  MPI_Comm               mpi_communicator(MPI_COMM_WORLD);
  ComputationalDomain<3> computational_domain(mpi_communicator);
  deal2lkit::ParameterAcceptor::initialize("test05_cad_import_and_global_refinement.prm","used.prm");
  computational_domain.read_domain();
  computational_domain.read_cad_files();
  computational_domain.assign_manifold_projectors(1.0e-6);


  if (refine==REFINE::GLOBAL)
    computational_domain.getTria().refine_global(3);
  else if (refine==REFINE::ASPECT_RATIO)
    computational_domain.aspect_ratio_refinement(4);


  std::string   filename0 = ("meshResult.inp");
  std::ofstream logfile0 (filename0.c_str ());
  GridOut       grid_out0;
  grid_out0.write_ucd (computational_domain.getTria(), logfile0);

}
