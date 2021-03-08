

#include "driver.h"
#include <iostream>

int
main(int argc, char *argv[])
{
  std::cout << "hello hello \n";


 try
   {
      for (int i=0;i<argc;++i)
        std::cout << "Input argument " << std::to_string(i) << ": " << argv[i] << "\n";



      unsigned int threads;
      if (argc == 1)
        threads = numbers::invalid_unsigned_int;
      else
        threads = atoi(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, threads);

      std::string pname =
        "parameters_bem_" + std::to_string(DEAL_II_DIMENSION) + ".prm";
      std::string pname2 =
        "used_parameters_bem_" + std::to_string(DEAL_II_DIMENSION) + ".prm";

      if (argc == 3)
        pname = argv[2];

      Driver<DEAL_II_DIMENSION> driver;
      deal2lkit::ParameterAcceptor::initialize(pname, pname2);

      driver.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
