#include <boost/program_options.hpp>
#include <deal.II/base/parameter_handler.h>

void exit_if_mpi (int argc, char *argv[]);

boost::program_options::variables_map
parse_cmdline_options (int argc, char *argv[]);

void parse_parameters_file (dealii::ParameterHandler &parameters,
                            std::string prm_file_name);
