#include "general.hpp"
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <boost/program_options.hpp>
#include <iostream>
#include <filesystem>

// using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;
namespace po = boost::program_options;

// Check if this app is runnint in parellel with "mpirun".
// Exit with error if true.
void exit_if_mpi (char** &argv) // name_of_executable
{
  if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) != 1)
  {
    if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      std::cerr << "This program can only be run in serial, use "
                << argv[0] << endl;
    }
    std::exit(1);
  }
}

// Parse command-line options.
po::variables_map parse_cmdline_options (int &argc, char **&argv)
{
  const unsigned int default_dim {2};
  const auto prm_file_name = "parameters.prm";

  // Allowd input options.
  po::options_description possible_options{"Options"};
  possible_options.add_options()
    ("dim,d", po::value<unsigned int>()->default_value(default_dim),
              "The dimension of the space.")
    ("parameters,p", po::value<std::string>()->default_value(prm_file_name),
                     "File with input parameters.");

  // Hidden options.
  po::options_description hidden_options{"Hidden options"};
  hidden_options.add_options()
    ("order,o", po::value<unsigned int>()->default_value(1),
                "Finite element order (polynomial degree).")
    ("mesh,m", po::value<std::string>(), "Input mesh file.")
    ("help,h", "Show help");

  po::options_description all_options{"All options"};
  all_options.add(possible_options)
             .add(hidden_options);

  // Positional options.
  po::positional_options_description positional_options;
  positional_options.add("parameters", -1);

  po::command_line_parser parser {argc, argv};
  parser.options(all_options)
        .positional(positional_options)
        .style( po::command_line_style::default_style |
                // Allow single dash for long options: -foo=10
                po::command_line_style::allow_long_disguise |
                po::command_line_style::allow_slash_for_short );
  po::parsed_options parsed_options = parser.run();

  po::variables_map options;
  po::store(parsed_options, options);
  po::notify(options);

  if(options.count("help")) {
    cout << possible_options << endl;
    std::exit(1);
  }

  cout << "Options used:" << endl
       << "  --parameters: " << options["parameters"].as<std::string>() << endl;
  // if (options.count("dim"))
  //   cout << "  --dim " << options["dim"].as<unsigned int>() << endl;
  // if (options.count("order"))
  //   cout << "  --order " << options["order"].as<unsigned int>() << endl;
  // if (options.count("mesh"))
  //   cout << "  --mesh: " << options["mesh"].as<std::string>() << endl;
  cout << endl;

  return options;
}

// Parse parameters file
void parse_parameters_file (dealii::ParameterHandler &parameters,
                            std::string prm_file_name)
{
  parameters.declare_entry(
      "Global mesh refinement steps",
      "5",
      dealii::Patterns::Integer(0, 20),
      "The number of times the 1-cell coarse mesh should be refined globally "
      "for our computations.");
  parameters.declare_entry(
      "Finite element order",
      "1",
      dealii::Patterns::Integer(1, 5),
      "Polynomial degree.");
  parameters.declare_entry(
      "Number of eigenvalues/eigenfunctions",
      "5",
      dealii::Patterns::Integer(0, 100),
      "The number of eigenvalues/eigenfunctions to be computed.");
  parameters.declare_entry(
      "Potential",
      "0",
      dealii::Patterns::Anything(),
      "A functional description of the potential.");
  parameters.parse_input( prm_file_name );
}


// vim: ts=2 sts=2 sw=2 fdm=marker
