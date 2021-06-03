#include "parameters_parsing.hpp"
#include "model.hpp"
#include "dft.hpp"

#include <deal.II/base/parameter_acceptor.h>
#include <boost/program_options.hpp>
#include <iostream>
#include <filesystem>

// using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;
namespace po = boost::program_options;


/// Build on top of the Boost.Program_options library.
po::variables_map parse_cmdline_options (int &argc, char** &argv)
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

  // Positional options. I.e. options that passed without any flags.
  // In this key it is only option 'parameters'.
  po::positional_options_description positional_options;
  positional_options.add("parameters", -1);

  po::command_line_parser parser {argc, argv};
  parser.options(all_options)
        .positional(positional_options)
        .style( po::command_line_style::default_style |
                po::command_line_style::allow_slash_for_short |
                po::command_line_style::allow_long_disguise );
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



template<int dim>
void parse_parameters (int &argc, char** &argv)
{
  // Parse command-line options.
  po::variables_map options = parse_cmdline_options (argc, argv);
  auto prm_file_name = options["parameters"].as<std::string>();

  Model<dim>::Parameters::get_parameters();
  DFT_Parameters::get_parameters();

  // Parse parameters file
  dealii::ParameterAcceptor::initialize(
      prm_file_name,
      "used_parameters.prm",
      dealii::ParameterHandler::Short |
        dealii::ParameterHandler::KeepDeclarationOrder,
      dealii::ParameterAcceptor::prm,
      dealii::ParameterHandler::KeepDeclarationOrder
  );
}



// {
//   parameters.declare_entry(
//       "Potential",
//       "0",
//       dealii::Patterns::Anything(),
//       "A functional description of the potential.");
//   parameters.parse_input( prm_file_name );
// }



/*--------------- Explicit templates instantiation --------------------------*/

template void parse_parameters<1>(int &argc, char** &argv);
template void parse_parameters<2>(int &argc, char** &argv);
template void parse_parameters<3>(int &argc, char** &argv);

// vim: ts=2 sts=2 sw=2
