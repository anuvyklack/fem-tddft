#include "general.hpp"
#include "model.hpp"
#include "ground_state.hpp"
// #include "time_dependent.hpp"

#include <deal.II/base/mpi.h>

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>

// using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  // const unsigned int dim {2};

  // Check if this app is runnint in parellel with "mpirun". Exit with error if true.
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  exit_if_mpi(argv);

  // Parse command-line options.
  po::variables_map options = parse_cmdline_options (argc, argv);
  // Parse parameters file
  dealii::ParameterHandler parameters;
  parse_parameters_file (parameters, options["parameters"].as<std::string>());



  cout << endl << "Job done." << endl;

  return 0;
}


// vim: ts=2 sts=2 sw=2 fdm=marker
