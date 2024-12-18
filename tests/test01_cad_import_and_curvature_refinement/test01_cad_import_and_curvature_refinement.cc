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
#include "computational_domain.h"

#include <filesystem>



void write_mesh_to_file(ComputationalDomain<3> &computational_domain, std::filesystem::path &filename);

int main(int argc, char **argv)
{
  std::string           this_file_name{__FILE__};
  std::filesystem::path cwd = std::filesystem::path(this_file_name).parent_path();
  std::filesystem::current_path(cwd);
  std::filesystem::path data_path = std::filesystem::path("../test_data/NACA_FILES/");
  std::filesystem::path output_path =
    cwd / std::filesystem::path("../output") / std::filesystem::path(this_file_name).stem();


  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  // initlog();
  MPI_Comm               mpi_communicator(MPI_COMM_WORLD);
  ComputationalDomain<3> computational_domain(mpi_communicator);

  deal2lkit::ParameterAcceptor::initialize(cwd / std::filesystem::path("parameters_bem_3.prm"), "used.prm");
  computational_domain.read_domain(data_path);
  computational_domain.read_cad_files(data_path);
  computational_domain.refine_and_resize(cwd);

  if (std::filesystem::exists(output_path) == false)
  {
    std::filesystem::create_directories(output_path);
  }

  std::filesystem::path filename = output_path / std::filesystem::path("meshResult.inp");

  write_mesh_to_file(computational_domain, filename);
}

void write_mesh_to_file(ComputationalDomain<3> &computational_domain, std::filesystem::path &filename)
{
  auto &mesh = computational_domain.getTria();
  {
    std::ofstream logfile(filename.c_str());
    logfile << std::setprecision(16);
    GridOut grid_out;
    grid_out.write_ucd(mesh, logfile);
  }
}
