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

#include <BRepBuilderAPI_MakeWire.hxx>


#include "./../tests.h"
#include "Body.h"
#include "WireUtil.h"
#include "computational_domain.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  // initlog();
  MPI_Comm               mpi_communicator(MPI_COMM_WORLD);
  ComputationalDomain<3> computational_domain(mpi_communicator);

  deal2lkit::ParameterAcceptor::initialize("test06_area_and_volume.prm", "used.prm");

  computational_domain.read_domain();

  Body             body;
  std::vector<int> wl = {11, 12};
  body.setWaterlineIndices(wl);

  auto shapes = computational_domain.cad_curves;

  BRepBuilderAPI_MakeWire wirebuilder;
  wirebuilder.Add(TopoDS::Wire(shapes[0]));
  wirebuilder.Add(TopoDS::Wire(shapes[1]));
  if (wirebuilder.IsDone())
  {
    double x0 = 100.0;
    double y0 = 0.0;
    double S0, Sx, Sy, Sxx, Sxy, Syy;
    WireUtil::surfaceMoments(wirebuilder.Wire(), x0, y0, S0, Sx, Sy, Sxx, Sxy, Syy);
    std::cout << "x0 : " << x0 << std::endl;
    std::cout << "y0 : " << y0 << std::endl;
    std::cout << "S0 : " << S0 << std::endl;
    std::cout << "Sx : " << Sx << std::endl;
    std::cout << "Sy : " << Sy << std::endl;
    std::cout << "Sxx : " << Sxx << std::endl;
    std::cout << "Sxy : " << Sxy << std::endl;
    std::cout << "Syy : " << Syy << std::endl;
  }
}
