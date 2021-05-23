#include "utilities.hpp"
#include "parameters_parsing.hpp"
#include "output_results.hpp"
// #include "output_results_2.hpp"
#include "model.hpp"
#include "dft.hpp"
#include "kohn_sham.hpp"
#include "hartree.hpp"
// #include "time_dependent.hpp"

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>

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
    DFT<dim> dft {model, dft_parameters};

    dft.run();


    // dealii::Vector<double> delta_density;
    // dealii::Vector<double> delta_hartree_potential;
    //
    // dft.solve_Kohn_Sham_problem();
    // dft.calculate_density();
    // dft.solve_Hartree_problem();
    //
    // delta_density = dft.data.density;
    // delta_hartree_potential = dft.data.hartree_potential;
    //
    // dft.solve_Kohn_Sham_problem();
    // dft.calculate_density();
    // dft.solve_Hartree_problem();
    //
    // delta_density -= dft.data.density;
    // delta_hartree_potential -= dft.data.hartree_potential;
    //
    // cout << delta_density.linfty_norm() << endl
    //      << delta_hartree_potential.linfty_norm() << endl;


    // Hartree<dim> hartree {model, dft.parameters, dft.data};
    // hartree.setup_system();
    // hartree.assemble_system();

    ResultsOutput<dim> data_out {model};
    data_out.add_data_vector(dft.data.kohn_sham_orbitals, "Kohn_Sham_orbitals");
    data_out.add_data_vector(dft.data.density, "density");
    data_out.add_data_vector(dft.data.hartree_potential, "hartree_potential");
    data_out.write("dft_results");

    // cout << "Calculate density." << endl;
    // dft.calculate_density();
    //
    // cout << "Solve Hartree problem." << endl;
    // dft.solve_Hartree_problem();

    // boost::regex_replace(str, boost::regex("[' ']{2,}"), " ");

    // model.output_ground_states();
    // model.save_to_file();  // Serialize ground state for futher calculations.

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
