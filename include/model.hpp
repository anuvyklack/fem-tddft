#ifndef MODEL_HPP
#define MODEL_HPP

#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/grid_out.h>
#include <fstream>
#include <filesystem>

enum ProblemType {
  GROUND_STATE,
  TIME_DEPENDENT
};

enum MeshType {
  SHELL,
  QUANTUM_WELL
};


// class DimensionlessModel;
class BaseProblem; // forward declaration

class DimensionlessModel
{
public:
  virtual void set_problem_type (ProblemType type) = 0;
  virtual void use_default_mesh (MeshType type) = 0;
  virtual void save_to_file (std::string file_name) const = 0;
  virtual void load_from_file() = 0;
  virtual void save_mesh() const = 0;

  // Getters and Setters ///////////////////////////////////////////////////////
  virtual void set_fe (std::unique_ptr <dealii::FiniteElement<1>> fe_to_set) { (void)fe_to_set; };
  virtual void set_fe (std::unique_ptr <dealii::FiniteElement<2>> fe_to_set) { (void)fe_to_set; };
  virtual void set_fe (std::unique_ptr <dealii::FiniteElement<3>> fe_to_set) { (void)fe_to_set; };
  virtual void get_fe(){};
  // virtual void = 0;
  // virtual void = 0;
  // virtual void = 0;
  // //////1/////////2/////////3/////////4/////////5/////////6/////////7////////

};


template <int dim>
class Model : public DimensionlessModel
// class Model
{
public:
  Model (dealii::ParameterHandler &parameters);
  void set_problem_type (ProblemType type);
  void use_default_mesh(MeshType type);

  void save_to_file ( std::string file_name = "data_for_restore" ) const;
  void load_from_file();
  void save_mesh() const;

  dealii::ParameterHandler &parameters;

  std::string name;
  std::filesystem::path results_path;

  dealii::Triangulation<dim> mesh;
  std::unique_ptr <dealii::FiniteElement<dim>> fe;
  dealii::DoFHandler<dim> dof_handler {mesh};

  std::vector<dealii::Vector<double>> ground_states;

  // Getters and Setters ///////////////////////////////////////////////////////
  void set_fe (std::unique_ptr <dealii::FiniteElement<dim>> fe_to_set) override;
  void get_fe();
  // //////1/////////2/////////3/////////4/////////5/////////6/////////7////////

// private:
  std::unique_ptr <BaseProblem> problem;
  // BaseProblem *problem;
};


class BaseProblem
{
public:
  virtual void run() = 0;
};

#endif
