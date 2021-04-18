#include "model.hpp"
#include "ground_state.hpp"
// #include "time_dependent.hpp"
#include <deal.II/grid/grid_generator.h>

using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;

template<int dim>
Model<dim>::Model (ParameterHandler &parameters)
  : parameters(parameters)
{}

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
      cout << "time dependent" << endl;
      // problem = std::make_unique<TDSE<dim>>(*this);
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
        const double inner_radius = 0.5,
                     outer_radius = 1.0;
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
void Model<dim>::save_to_file ( std::string file_name ) const
{
  fs::path file_path = results_path / file_name;
  std::ofstream file(file_path);
  // boost::archive::text_oarchive archive(file);
  boost::archive::binary_oarchive archive(file);
  archive << fe->degree << mesh << dof_handler << ground_states;
}

template <int dim>
void Model<dim>::load_from_file()
{

}

template <int dim>
void Model<dim>::save_mesh() const
{
  fs::path file_path = results_path / (name + ".vtk");
  std::ofstream file(file_path);
  GridOut grid_out;
  grid_out.write_vtk(mesh, file);
}

// Getters and Setters ///////////////////////////////////////////////////////
template <int dim>
void Model<dim>::set_fe (std::unique_ptr<FiniteElement<dim>> fe_to_set)
{
  fe = std::move(fe_to_set);
  cout << fe->get_name() << endl;
}

template <int dim>
void Model<dim>::get_fe() {}
// //////1/////////2/////////3/////////4/////////5/////////6/////////7////////

// Explicit template instantiation
template class Model<1>;
template class Model<2>;
template class Model<3>;


// vim: ts=2 sts=2 sw=2 fdm=marker
