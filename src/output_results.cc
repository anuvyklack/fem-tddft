#include "output_results.hpp"
#include "models/base_model.hpp"

#include <deal.II/numerics/data_out.h>
#include <magic_enum.hpp>
#include <iostream>
#include <filesystem>

// using namespace dealii;
// using std::cin, std::cout, std::endl;
// namespace fs = std::filesystem;
//
//
// template <int dim>
// ResultsOutput<dim>::ResultsOutput (const Model<dim> & model)
//   : model(model)
// {
//   data_out.attach_dof_handler(model.dof_handler);
//   set_output_format("vtu");
// }
//
//
//
// template <int dim>
// void ResultsOutput<dim>::set_output_format (std::string format_name)
// {
//   auto format =
//       magic_enum::enum_cast<dealii::DataOutBase::OutputFormat>(format_name);
//   if (format.has_value())
//     output_format = format.value();
//   else
//     AssertThrow(false, dealii::ExcMessage("Invalide output format!"));
//
//   // Filename extension.
//   extension = "." + format_name;
// }
//
//
//
// template <int dim>
// void ResultsOutput<dim>::add_fe_function (dealii::Vector<double> function,
//                                           std::string function_name)
// {
//   data_out.add_data_vector(function, function_name);
//   data_out.build_patches();
// }
//
//
//
// template <int dim>
// template <typename vec_type>
// void ResultsOutput<dim>::add_fe_function(
//     std::vector<dealii::Vector<vec_type>> functions,
//     std::string functions_name
// )
// {
//   for (unsigned int i = 0; i < functions.size(); ++i)
//     data_out.add_data_vector(
//         functions[i],
//         functions_name + Utilities::int_to_string(i)
//     );
//   data_out.build_patches();
// }
//
//
//
// template <int dim>
// void ResultsOutput<dim>::write (std::string file_name)
// {
//   fs::path file_path = model.results_path / (file_name + extension);
//   std::ofstream output(file_path);
//   data_out.write(output, output_format);
//
// }
//
//
// /*--------------- Explicit templates instantiation --------------------------*/
//
// template<> template<>
// void ResultsOutput<1>::add_fe_function(
//     std::vector<dealii::Vector<double>> functions,
//     std::string functions_name
// );
//
// template<> template<>
// void ResultsOutput<2>::add_fe_function(
//     std::vector<dealii::Vector<double>> functions,
//     std::string functions_name
// );
//
// template<> template<>
// void ResultsOutput<3>::add_fe_function(
//     std::vector<dealii::Vector<double>> functions,
//     std::string functions_name
// );
//
// template class ResultsOutput<1>;
// template class ResultsOutput<2>;
// template class ResultsOutput<3>;

// vim: ts=2 sts=2 sw=2
