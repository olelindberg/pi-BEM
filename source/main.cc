

#include <boost/program_options.hpp>

#include <exception>
#include <iostream>

#include "driver.h"

#include <filesystem>

class Options
{
public:
  Options(int argc, char *argv[])
  {
    _desc = std::make_shared<boost::program_options::options_description>("pi-BEM options");
    _desc->add_options()("help", "produce help message");
    _desc->add_options()("input-path", boost::program_options::value<std::string>(), "Input path");
    _desc->add_options()("threads", boost::program_options::value<int>(), "Number of threads");

    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, *_desc),
                                  _vm);
    boost::program_options::notify(_vm);

    if (_vm.count("help"))
      std::cout << *_desc << "\n";

    if (_vm.count("input-path"))
      _input_path = _vm["input-path"].as<std::string>();
  }

  std::string get_input_path()
  {
    return _input_path;
  };

private:
  std::shared_ptr<boost::program_options::options_description> _desc;
  boost::program_options::variables_map                        _vm;
  std::string                                                  _input_path;
};

int main(int argc, char *argv[])
{
  try
  {
    Options options(argc, argv);

    unsigned int threads;
    if (argc == 1)
      threads = numbers::invalid_unsigned_int;
    else
      threads = atoi(argv[1]);
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, threads);
    std::cout << "MPI initialized with " << threads << " threads" << std::endl;

    auto output_path = std::filesystem::path(options.get_input_path()).append("output").string();
    if (!std::filesystem::exists(output_path))
    {
      std::cout << "Creating output directory: " << output_path << "\n";
      std::filesystem::create_directory(output_path);
    }
    std::string pname = std::filesystem::path(options.get_input_path())
                          .append("parameters_bem_" + std::to_string(DEAL_II_DIMENSION) + ".prm")
                          .string();
    std::string pname2 =
      std::filesystem::path(output_path)
        .append("used_parameters_bem_" + std::to_string(DEAL_II_DIMENSION) + ".prm")
        .string();

    if (argc == 3)
      pname = argv[2];

    Driver<DEAL_II_DIMENSION> driver;
    deal2lkit::ParameterAcceptor::initialize(pname, pname2);

    driver.run(options.get_input_path(), output_path);
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------" << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------" << std::endl;

    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------" << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------" << std::endl;
    return 1;
  }

  return 0;
}
