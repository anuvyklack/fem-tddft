#ifndef MODEL_HPP
#define MODEL_HPP

// #include "ground_state.hpp"

#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/grid_out.h>
#include <fstream>
#include <filesystem>

/// @brief The problem types that can be attached to a Model to solve.
enum ProblemType {
  GROUND_STATE,
  HARTREE,
  TIME_DEPENDENT
};

/// @brief The built-in meshes types, that can be used instead of specifying an input
/// file with mesh to read.
enum MeshType {
  SHELL,
  QUANTUM_WELL
};

// forward declaration
class BaseProblem;
// template<int dim> class EigenvalueProblem;
// template<int dim> class TDSE;

/**
 * @brief The main class that represents the considering finite-element model.
 */
template <int dim>
class Model
{
public:
  Model (const dealii::ParameterHandler &parameters);
  void set_problem_type (ProblemType type);
  void use_default_mesh(MeshType type);

  /// @brief Save the current Model state into file.
  void save_to_file ( std::string file_name = "data_for_restore" ) const;

  /// @brief Restore the saved Model state from file.
  void load_from_file ( std::string file_name = "data_for_restore" );

  void save_mesh() const;
  void output_ground_states() const;

  const dealii::ParameterHandler &parameters;

  std::string name;
  std::filesystem::path results_path;

  dealii::Triangulation<dim> mesh;
  dealii::DoFHandler<dim> dof_handler {mesh};

  std::vector<dealii::Vector<double>> stationary_states; // стационарные состояния

  std::unique_ptr <BaseProblem> problem;

  // Getters and Setters //////////////////////////////////////////////////////
  const dealii::FiniteElement<dim>& get_fe() const { return *fe; }
  // ////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

  // friend class EigenvalueProblem<dim>;
  // friend class TDSE<dim>;

private:
  void set_fe (const unsigned int &order);

  std::unique_ptr <dealii::FiniteElement<dim>> fe;
};


class BaseProblem
{
public:
  virtual void run() = 0;
};

#endif
// vim: ts=2 sts=2 sw=2
