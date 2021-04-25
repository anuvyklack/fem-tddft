#ifndef GENERAL_HPP
#define GENERAL_HPP

#include <boost/program_options.hpp>
#include <deal.II/base/parameter_handler.h>

/// @brief Check if this app is running in parellel with "mpirun".
/// Exit the application with error if true.
void exit_if_mpi (char** &argv);

/// @brief Parse command-line options.
boost::program_options::variables_map
  parse_cmdline_options (int &argc, char** &argv);

/// @brief Parse parameters file.
void parse_parameters_file (const std::string &prm_file_name,
                            dealii::ParameterHandler &parameters);

#endif
// vim: ts=2 sts=2 sw=2
