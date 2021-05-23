#ifndef OUTPUT_RESULTS_HEADER
#define OUTPUT_RESULTS_HEADER

#include "model.hpp"
#include <deal.II/numerics/data_out.h>
#include <magic_enum.hpp>
#include <iostream>
#include <filesystem>


template<int dim>
class ResultsOutput
{
public:
  ResultsOutput (const Model<dim> & model);
  void set_output_format (std::string format_name);
  void add_data_vector (const dealii::Vector<double> &function,
                        std::string function_name);

  template <typename VectorType>
  void add_data_vector (const std::vector<VectorType> &functions,
                        std::string function_name);

  /// Write to file
  void write (std::string file_name);

private:
  dealii::DataOut<dim> data_out;

  const Model<dim> & model;

  dealii::DataOutBase::OutputFormat output_format;
  std::string extension; ///< Filename extension;

};


/*--------------- Tlmplate and inline functions -----------------------------*/

template <int dim>
ResultsOutput<dim>::ResultsOutput (const Model<dim> & model)
  : model(model)
{
  data_out.attach_dof_handler(model.dof_handler);
  set_output_format("vtu");
}



template <int dim>
void ResultsOutput<dim>::set_output_format (std::string format_name)
{
  auto format =
      magic_enum::enum_cast<dealii::DataOutBase::OutputFormat>(format_name);
  if (format.has_value())
    output_format = format.value();
  else
    AssertThrow(false, dealii::ExcMessage("Invalide output format!"));

  // Filename extension.
  extension = "." + format_name;
}



template <int dim>
void ResultsOutput<dim>::add_data_vector(
    const dealii::Vector<double> &function,
    std::string function_name
)
{
  data_out.add_data_vector(function, function_name);
  // data_out.build_patches();

  // Vector<double> projected_potential(dof_handler.n_dofs());
  // {
  //   FunctionParser<dim> potential;
  //   potential.initialize(FunctionParser<dim>::default_variable_names(),
  //                        parameters.get("Potential"),
  //                        typename FunctionParser<dim>::ConstMap());
  //   VectorTools::interpolate(dof_handler, potential, projected_potential);
  // }
  // data_out.add_data_vector(projected_potential, "interpolated_potential");

}



template <int dim>
template <typename VectorType>
void ResultsOutput<dim>::add_data_vector(
    const std::vector<VectorType> &functions,
    std::string functions_name
)
{
  for (unsigned int i = 0; i < functions.size(); ++i)
    data_out.add_data_vector(
        functions[i],
        functions_name + dealii::Utilities::int_to_string(i)
    );

  // data_out.build_patches();

}



template <int dim>
void ResultsOutput<dim>::write (std::string file_name)
{
  data_out.build_patches();

  std::filesystem::path file_path = model.results_path / (file_name + extension);
  std::ofstream output(file_path);
  data_out.write(output, output_format);
}



#endif // OUTPUT_RESULTS_HEADER
// vim: ts=2 sts=2 sw=2
