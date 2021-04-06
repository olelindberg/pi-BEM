#include <iostream>
#include "driver.h"

int main(int argc, char *argv[])
{
  unsigned int                     threads = numbers::invalid_unsigned_int;
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, threads);
  std::string                      pname  = "parameters_bem_3.prm";
  std::string                      pname2 = "used_parameters_bem_3.prm";
  MPI_Comm                         mpi_communicator(MPI_COMM_WORLD);
  ComputationalDomain<3>           computational_domain(mpi_communicator);
  deal2lkit::ParameterAcceptor::initialize(pname, pname2);
  std::cout << "main00" << std::endl;
  computational_domain.read_domain();  
  std::cout << "main01" << std::endl;
  bool with_double_nodes = false;
  computational_domain.make_edges_conformal(with_double_nodes);
  std::cout << "main02" << std::endl;
  return 0;
}
