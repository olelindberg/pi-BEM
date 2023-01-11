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

#include <deal.II/grid/grid_tools.h>

#include "./../tests.h"
#include "computational_domain.h"

enum class REFINE
{
  GLOBAL,
  ASPECT_RATIO
};

int
main(int argc, char **argv)
{
  auto refine = REFINE::GLOBAL;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  // initlog();
  MPI_Comm               mpi_communicator(MPI_COMM_WORLD);
  ComputationalDomain<3> computational_domain(mpi_communicator);
  deal2lkit::ParameterAcceptor::initialize("test15_cad_bathymetry.prm", "used.prm");
  computational_domain.read_domain();

  // --------------------------------------------------------------------------
  // Rotate and translate:
  // --------------------------------------------------------------------------
  // double pi    = std::acos(-1.0);
  // double angle = 0.9 * pi;
  // GridTools::rotate(angle, 2, computational_domain.tria);

  // Tensor<1, 3> translation;
  // translation[0] = 582883.0;
  // translation[1] = 6314420.0;
  // translation[2] = 0.0;
  // GridTools::shift(translation, computational_domain.tria);

  // --------------------------------------------------------------------------
  // Project to manifold:
  // --------------------------------------------------------------------------
  // auto &tria = computational_domain.tria;
  // for (auto cell = tria.begin_active(); cell != tria.end(); ++cell)
  // {
  //   std::vector<Point<3>> vertices(GeometryInfo<2>::vertices_per_cell);
  //   for (unsigned int i = 0; i < GeometryInfo<2>::vertices_per_cell; ++i)
  //     vertices[i] = cell->vertex(i);

  //   for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
  //     cell->vertex(v) = cell->get_manifold().project_to_manifold(vertices, cell->vertex(v));
  // }

  computational_domain.getTria().refine_global(2);

  std::string   filename0 = ("meshResult.inp");
  std::ofstream logfile0(filename0.c_str());
  logfile0 << std::setprecision(16);
  GridOut grid_out0;
  grid_out0.write_ucd(computational_domain.getTria(), logfile0);
}
