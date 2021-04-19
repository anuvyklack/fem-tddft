#include "model.hpp"
#include "ground_state.hpp"
#include "time_dependent.hpp"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/numerics/data_out.h>

using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;

template<int dim>
Model<dim>::Model (ParameterHandler &parameters)
  : parameters(parameters)
{
  set_fe( parameters.get_integer("Finite element order") );
}

template<int dim>
void Model<dim>::set_fe (const unsigned int &order)
{
  fe = std::make_unique <dealii::FE_Q<dim>> (order);
}

template<int dim>
void Model<dim>::set_problem_type (ProblemType type)
{
  switch (type)
  {
    case GROUND_STATE:
      // cout << "ground state" << endl;
      problem = std::make_unique<EigenvalueProblem<dim>>(*this);
      break;
    case TIME_DEPENDENT:
      // cout << "time dependent" << endl;
      problem = std::make_unique<TDSE<dim>>(*this);
      break;
  };
}

template <int dim>
void Model<dim>::use_default_mesh (MeshType type)
{
  switch (type)
  {
    case SHELL:
    {
      name = "shell";
      const double inner_radius = 0.5, outer_radius = 1.0;

      // Empty Point<dim>{} constructor creates a point at the origin.
      GridGenerator::hyper_shell(mesh, Point<dim>{},
                                 inner_radius, outer_radius,
                                 24);  // 24, 48, 96, 192*2^m
      break;
    }
    case QUANTUM_WELL:
      name = "quantum_well";
      GridGenerator::hyper_cube(mesh, -1, 1);
      break;
  }
}

template <int dim>
void Model<dim>::save_to_file (std::string file_name) const
{
  fs::path file_path = results_path / file_name;
  std::ofstream file(file_path);
  // boost::archive::text_oarchive archive(file);
  boost::archive::binary_oarchive archive(file);
  archive << fe->degree << mesh << dof_handler << ground_states;
}

template <int dim>
void Model<dim>::load_from_file (std::string file_name)
{
  fs::path file_path = results_path / file_name;
  std::ifstream file(file_path);
  boost::archive::binary_iarchive archive(file);
  unsigned int fe_degree;
  archive >> fe_degree;
  set_fe(fe_degree);
  // dof_handler.distribute_dofs(fe);
  archive >> mesh >> dof_handler >> ground_states;
}

template <int dim>
void Model<dim>::save_mesh() const
{
  fs::path file_path = results_path / (name + ".vtk");
  std::ofstream file(file_path);
  GridOut grid_out;
  grid_out.write_vtk(mesh, file);
}

template <int dim>
void Model<dim>::output_ground_states() const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  for (unsigned int i = 0; i < ground_states.size(); ++i)
    data_out.add_data_vector(ground_states[i],
                             std::string("grond_state_") +
                                 Utilities::int_to_string(i));
  data_out.build_patches();

  fs::path file_path = results_path / "ground_states.vtu";
  std::ofstream output(file_path);
  data_out.write_vtu(output);
}


// Explicit template instantiation
template class Model<1>;
template class Model<2>;
template class Model<3>;


// template <int dim>
// void print_mesh_info(const Triangulation<dim> &mesh)
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
// }


// vim: ts=2 sts=2 sw=2 fdm=marker
