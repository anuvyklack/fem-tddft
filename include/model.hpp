#ifndef MODEL_HEADER
#define MODEL_HEADER

#include "parameters_parsing.hpp"
#include "double_output_stream.hpp"
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/timer.h>
#include <boost/any.hpp>
#include <fstream>
#include <filesystem>


// /// @brief The problem types that can be attached to a Model to solve.
// enum class ProblemType
// {
//   Stationary_Schrodinger, ///< Stationary Schrodinger equation.
//   Hartree,
//   DFT,
//   DFT_Uniform, ///< Uniform mesh refinment.
//   DFT_Adaptive, ///< Adaptive mesh refinment.
//   Time_Dependent ///< Time dependent Schrodinger equation.
// };

/// @brief The built-in meshes types.
/// Can be used instead of specifying an input file with mesh to read.
enum BuiltinMesh
{
  quantum_well,
  shell,
};

/**
 * @brief The main class that represents the considering finite-element model.
 */
template <int dim>
class Model
{
public:
  /// The Meyers Singleton.
  class Parameters : public dealii::ParameterAcceptor
  {
  public:
    static Parameters & get_parameters()
    {
      static Model<dim>::Parameters prm;
      return prm;
    }

    /**
     * If @p true then @link which_mesh_to_use @endlink attribute stores
     *  value of type * @link BuiltinMesh @endlink.
     * If @p false @link which_mesh_to_use @endlink stores path to mesh file.
     */
    bool use_built_in_mesh;

    /**
     * Stores either mesh file path or one of built-in meshes identifiers
     * (value of type @link BuiltinMesh @endlink). Which one is stored
     * depends on @link use_built_in_mesh @endlink property.
     */
    boost::any mesh_to_use;

    unsigned int global_mesh_refinement_steps = 5;

    /// The finite element name to create finite element object from.
    std::string  fe_name;
    std::string  fe_type = "FE_Q";
    unsigned int fe_order = 1;

    std::string results_folder = "auto";

    unsigned int verbosity_level = 10;

    bool initialized = false;

  private:
    Parameters(); /* Constructor should be private! */
    Parameters (const Parameters& other) = delete;
    Parameters& operator=(const Parameters&) = delete;

    virtual void parse_parameters (dealii::ParameterHandler & prm) override;
  };

  Model (Parameters & parameters);

  ~Model();

  /// @brief Save the current Model state into file.
  void save_to_file (std::string file_name = "data_for_restore") const;

  /// @brief Restore the saved Model state from file.
  void load_from_file (std::string file_name = "data_for_restore");

  void output_mesh() const;
  // void output_stationary_states() const;

  /// @brief Returns the constant reference to @p fe attribute.
  inline const dealii::FiniteElement<dim>& get_fe() const { return *fe_ptr; }

  const Parameters & parameters;

  std::string name;
  std::filesystem::path results_path;

  double_ostream out;

  dealii::Triangulation<dim> mesh;
  dealii::DoFHandler<dim> dof_handler {mesh};

private:
  void set_mesh ();

  std::unique_ptr <dealii::FiniteElement<dim>> fe_ptr;

  std::ofstream output_file;
  std::ofstream log_file;
  // std::ofstream log_file {results_path / "log"};
};



class BaseProblem
{
protected:
  // dealii::TimerOutput monitor;
};


#endif // MODEL_HEADER
// vim: ts=2 sts=2 sw=2
