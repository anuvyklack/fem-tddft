#include "general.hpp"
#include "model.hpp"
// #include "stationary_schrodinger.hpp"
// #include "time_dependent.hpp"

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>

// using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;
namespace po = boost::program_options;

// int main(int argc, char *argv[])
int main(int argc, char** argv)
{
  const unsigned int dim {1};

  // Check if this app is runnint in parellel with "mpirun". Exit with error if true.
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  exit_if_mpi(argv);

  // Parse command-line options.
  po::variables_map options = parse_cmdline_options (argc, argv);
  // Parse parameters file
  dealii::ParameterHandler parameters;
  parse_parameters_file (options["parameters"].as<std::string>(), parameters);

  parameters.print_parameters(cout,
      dealii::ParameterHandler::PRM |
      // dealii::ParameterHandler::LaTeX |
      dealii::ParameterHandler::KeepDeclarationOrder
  );


  Model<dim> model {parameters};

  model.use_default_mesh(QUANTUM_WELL);
  // model.use_default_mesh(SHELL);

  model.results_path = fs::current_path() /
                       (model.name + "_" + std::to_string(dim) + "D_result");
  fs::create_directory(model.results_path);

  // cout << model.get_fe().get_name() << endl;

  model.set_problem_type(GROUND_STATE);
  model.problem->run();
  model.output_ground_states();
  model.save_to_file();  // Serialize ground state for futher calculations.

  model.set_problem_type(HARTREE);
  model.problem->run();

  // model.set_problem_type(TIME_DEPENDENT);
  // model.problem->run();

  // model.save_to_file();

  cout << endl << "Job done." << endl;

  return 0;
}


// vim: ts=2 sts=2 sw=2 fdm=marker
