#include "utilities.hpp"
#include <deal.II/base/mpi.h>
// #include <deal.II/base/parameter_handler.h>
// #include <boost/program_options.hpp>
#include <deal.II/numerics/data_out.h>
#include <iostream>
#include <filesystem>

using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;


void exit_if_mpi (char** &argv)  // std::string name_of_executable
{
  if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) != 1)
    {
      if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          std::cerr << "This program can only be run in serial, use "
                    << "./" << argv[0] << endl;
        }
      std::exit(1);
    }
}


/*------------------ Explicit templates instantiation -------------------*/

// vim: ts=2 sts=2 sw=2
