#ifndef PARAMETERS_PARSING_HEADER
#define PARAMETERS_PARSING_HEADER

#include <boost/program_options.hpp>
#include <deal.II/base/parameter_handler.h>

/// @brief Parse command-line options.
boost::program_options::variables_map
parse_cmdline_options (int &argc, char** &argv);

/// @brief Parse parameters file.
template<int dim>
void parse_parameters (int &argc, char** &argv);

#endif // PARAMETERS_PARSING_HEADER
// vim: ts=2 sts=2 sw=2
