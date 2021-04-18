#include "general.hpp"
#include "model.hpp"
#include "ground_state.hpp"
// #include "time_dependent.hpp"

#include <deal.II/base/mpi.h>
// #include <deal.II/grid/grid_generator.h>

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>

// using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  const unsigned int dim {2};

  // Check if this app is runnint in parellel with "mpirun".
  // Exit with error if true.
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  exit_if_mpi(argc, argv);

  // Parse command-line options.
  po::variables_map options = parse_cmdline_options (argc, argv);

  // Parse parameters file
  dealii::ParameterHandler parameters;
  parse_parameters_file (parameters, options["parameters"].as<std::string>());


  std::unique_ptr<DimensionlessModel> model;

  switch (dim)
  {
    case 1:
    {
      const unsigned int dim {1};
      break;
    }

    case 2:
    {
      const unsigned int dim {2};
      model = std::make_unique<Model<dim>> (parameters);
      model->set_fe( std::make_unique <dealii::FE_Q<dim>>(
                         parameters.get_integer("Finite element order")) );
      break;
    }

    case 3:
    {
      break;
    }
  }

  // Model<dim> model {parameters};
  //
  // model.fe = std::make_unique <dealii::FE_Q<dim>> (
  //                parameters.get_integer("Finite element order"));
  //
  // model.use_default_mesh(QUANTUM_WELL);
  // // model.use_default_mesh(SHELL);
  //
  //
  // model.results_path = fs::current_path() /
  //                        (model.name + "_" + std::to_string(dim) + "D_result");
  // fs::create_directory(model.results_path);
  //
  // model.set_problem_type(GROUND_STATE);
  // // model.set_problem_type(TIME_DEPENDENT);
  //
  // // model.problem = std::make_unique<EigenvalueProblem<dim>>(model);
  // model.problem->run();
  //
  // // {
  // //   EigenvalueProblem<dim> ground_state(model);
  // //   ground_state.run();
  // // }
  //
  // // cout << model.ground_states.size() << endl;
  //
  // // model.save_to_file();

  cout << endl << "Job done." << endl;

  return 0;
}


// template <int dim>
// void print_mesh_info(const Triangulation<dim> &mesh) //{{{
// {
//   cout << "Mesh info:" << endl
//        << "  dimension: " << dim << endl
//        << "  number of cells: " << mesh.n_active_cells() << endl;
//   {
//     std::map<types::boundary_id, unsigned int> boundary_count;
//     for (const auto &face : mesh.active_face_iterators())
//       if (face->at_boundary())
//         boundary_count[face->boundary_id()]++;
//
//     cout << "  boundary indicators: ";
//     for (const std::pair<const types::boundary_id, unsigned int> &pair : boundary_count)
//       {
//         cout << pair.first << " (" << pair.second << " times) ";
//       }
//     cout << endl;
//   }
// } //}}}


// vim: ts=2 sts=2 sw=2 fdm=marker
