#ifndef GENERAL_HPP
#define GENERAL_HPP

#include <boost/program_options.hpp>
#include <deal.II/base/parameter_handler.h>

void exit_if_mpi (char **&argv); // std::string name_of_executable

boost::program_options::variables_map
parse_cmdline_options (int &argc, char **&argv);

void parse_parameters_file (dealii::ParameterHandler &parameters,
                            std::string prm_file_name);

#endif
// vim: ts=2 sts=2 sw=2
