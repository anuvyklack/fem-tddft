#include "utilities.hpp"
#include "parameters_parsing.hpp"
#include "output_results.hpp"
// #include "output_results_2.hpp"
#include "model.hpp"
#include "dft.hpp"
#include "kohn_sham.hpp"
#include "hartree.hpp"
// #include "time_dependent.hpp"
#include "external_potentials.hpp"

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/numerics/vector_tools.h>

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>

// using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;
namespace po = boost::program_options;

// template<int dim>
// typename Model<dim>::Parameters model_parameters;

// int main(int argc, char *argv[])
int main(int argc, char** argv)
{
try
  {
    const unsigned int dim {1};

    // Check if this app is runnint in parellel with "mpirun".
    // Exit with error if true.
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    exit_if_mpi(argv);

    parse_parameters<dim>(argc, argv);

    Model<dim>::Parameters & model_parameters =
                             Model<dim>::Parameters::get_parameters();
    Model<dim> model {model_parameters};

    DFT_Parameters & dft_parameters = DFT_Parameters::get_parameters();
    // PointCharge<dim> potential {1, dealii::Point<dim>()};

    // DFT<dim> dft {model, dft_parameters, potential};
    DFT<dim> dft {model, dft_parameters};
    dft.run();


    ResultsOutput<dim> data_out {model};
    data_out.add_data_vector(dft.kohn_sham_orbitals.wavefunctions, "Kohn_Sham_orbital");
    data_out.add_data_vector(dft.density, "density");
    data_out.add_data_vector(dft.hartree_potential, "hartree_potential");

    // dealii::Vector<double> fe_potential (model.dof_handler.n_dofs());
    // dealii::VectorTools::interpolate( model.dof_handler, potential, fe_potential );
    // data_out.add_data_vector(fe_potential, "point_charge");

    data_out.write("dft_results");

    // boost::regex_replace(str, boost::regex("[' ']{2,}"), " ");
  }
catch (std::exception &exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
catch (...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  cout << endl << "Job done." << endl;
  return 0;
}


// vim: ts=2 sts=2 sw=2 fdm=marker
