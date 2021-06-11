#include "models/base_model.hpp"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/data_out.h>
#include <magic_enum.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>

using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;



template <int dim>
Model<dim>::Parameters::Parameters()
  : ParameterAcceptor("/General")
{
  add_parameter(
      "Results folder",
      results_folder,
      "The folder to store results to. "
      "If set to 'auto' the results folder name will be created atomaticaly "
      "from mesh name.");

  // // declare Mesh parameter
  // enter_my_subsection(prm);
  // {
  //   // TODO: with C++20 this string can be 'constexpr'.
  //   std::string description =
  //       std::string { "Mesh to use. Should be path to mesh file or "
  //                     "one of built-in meshes names.\n"
  //                     "Built-in meshes:\n" }
  //       +
  //       []()->std::string
  //       {
  //         std::string result;
  //         // Cycle through 'std::array' with enum values
  //         // converted to 'std::view_string'.
  //         for (const auto& mesh_type : magic_enum::enum_names<BuiltinMesh>())
  //           result += "\t- " + static_cast<std::string>(mesh_type) + "\n";
  //         return result;
  //       }();
  //
  //   prm.declare_entry(
  //       "Mesh",
  //       "quantum_well",
  //       dealii::Patterns::Anything(),
  //       description);
  // }
  // leave_my_subsection(prm);

  add_parameter(
      "Finite element space",
      fe_type,
      "The finite element space to use. "
      "This option accepts any standard finite element name "
      "defined in the DEAL II library.");
  add_parameter(
      "Finite element order",
      fe_order,
      "Polynomial degree.");
  add_parameter(
      "Global mesh refinement steps",
      global_mesh_refinement_steps,
      "The number of times the 1-cell coarse mesh should be"
      "refined globally for our computations.");
  add_parameter(
      "Verbosity level",
      verbosity_level,
      "Set the verbosity level in the log file.");

  // parse_parameters_call_back.connect(
  //     std::bind( &Model<dim>::Parameters::set_mesh_type, this )
  // );

  parse_parameters_call_back.connect( [&](){initialized = true;} );
}



template <int dim>
void Model<dim>::Parameters::parse_parameters (ParameterHandler & prm)
{
  (void)prm;

  // // Parse Mesh parameter.
  // std::string mesh = prm.get("Mesh");
  //
  // // Check if string `mesh` is in array `enum_names<>()`.
  // if (boost::algorithm::any_of_equal(
  //         magic_enum::enum_names<BuiltinMesh>(), mesh ))
  //   {
  //     use_built_in_mesh = true;
  //     mesh_to_use = magic_enum::enum_cast<BuiltinMesh>( mesh ).value();
  //   }
  // else
  //   {
  //     use_built_in_mesh = false;
  //     // TODO: parse `prm.get("Mesh")` result as filepath
  //     // and store to `which_mesh_to_use`.
  //     AssertThrow(false, ExcMessage(
  //                 "This functionality is not implemented yet."));
  //   }

  // Parse finite element name.
  fe_name = fe_type + '(' + std::to_string(fe_order) + ')';

  // // Parse results folder name.
  // if (results_folder == "auto" && use_built_in_mesh)
  //   {
  //     results_folder = mesh + "_" + std::to_string(dim) + "D";
  //   }
  // else if (results_folder == "auto")
  //   {
  //     // TODO: get file name from the mesh file path.
  //     AssertThrow(false, ExcMessage(
  //                 "This functionality is not implemented yet."));
  //   }

  // Parse results folder name.
  if (results_folder != "auto")
    {
      results_folder = results_folder + "_" + std::to_string(dim) + "D";
    }
  else if (results_folder == "auto")
    {
      // TODO: get file name from the mesh file path.
      AssertThrow(false, ExcMessage(
                  "This functionality is not implemented yet."));
    }
}




template<int dim>
Model<dim>::Model (Parameters & parameters)
  : parameters(parameters),
    out(cout, output_file)
{
  AssertThrow(parameters.initialized, ExcNotInitialized())

  // Create the results folder.
  results_path = fs::current_path() / parameters.results_folder;
  fs::create_directory(results_path);

  // Open output file stream.
  output_file.open(results_path / "output.txt");

  // set_mesh();
  fe_ptr = FETools::get_fe_by_name<dim>( parameters.fe_name );

  // out << "Global mesh refinement steps: "
  //     << std::to_string( parameters.global_mesh_refinement_steps ) << endl
  //     << endl;
  // mesh.refine_global( parameters.global_mesh_refinement_steps );
  //
  // dof_handler.distribute_dofs(*fe_ptr);

  // Setup logging.
  log_file.open(results_path / "log");
  deallog.attach(log_file);
  deallog.depth_file( parameters.verbosity_level );
}



template<int dim>
Model<dim>::~Model()
{
  output_file.close();
  log_file.close();
}



// template <int dim>
// void Model<dim>::set_mesh ()
// {
//   if (parameters.use_built_in_mesh)
//     {
//       switch (boost::any_cast<BuiltinMesh>(parameters.mesh_to_use))
//         {
//           case shell:
//           {
//             const double inner_radius = 0.5, outer_radius = 1.0;
//
//             // Empty Point<dim>{} constructor creates a point at the origin.
//             GridGenerator::hyper_shell(mesh, Point<dim>{},
//                                        inner_radius, outer_radius,
//                                        24);  // 24, 48, 96, 192*2^m
//             break;
//           }
//           case quantum_well:
//           {
//             GridGenerator::hyper_cube(mesh, -1, 1);
//             break;
//           }
//         }
//     }
//   else
//     {
//       // TODO: load mesh form file.
//       AssertThrow(false, ExcMessage(
//                   "This functionality is not implemented yet."));
//     }
// }



template <int dim>
void Model<dim>::save_to_file (std::string file_name) const
{
  std::ofstream file (results_path / file_name);
  // boost::archive::text_oarchive archive(file);
  boost::archive::binary_oarchive archive (file);
  archive << fe_ptr->get_name() << mesh << dof_handler;
}



template <int dim>
void Model<dim>::load_from_file (std::string file_name)
{
  std::ifstream file (results_path / file_name);
  boost::archive::binary_iarchive archive (file);

  // Restore finite element object.
  std::string fe_name;
  archive >> fe_name;
  fe_ptr = FETools::get_fe_by_name<dim>(fe_name);

  // dof_handler.distribute_dofs(fe);
  archive >> mesh >> dof_handler;
}



template <int dim>
void Model<dim>::output_mesh() const
{
  GridOut::OutputFormat output_format = GridOut::vtu;
  // filename extension
  std::string ext = "." + static_cast<std::string>(
      magic_enum::enum_name<GridOut::OutputFormat>(output_format) );

  fs::path file_path = results_path / ("mesh" + ext);
  std::ofstream file(file_path);
  GridOut grid_out;
  grid_out.write(mesh, file, output_format);
}



// template <int dim>
// void Model<dim>::output_stationary_states() const
// {
//   DataOut<dim> data_out;
//   data_out.attach_dof_handler(dof_handler);
//
//   for (unsigned int i = 0; i < stationary_states.size(); ++i)
//     data_out.add_data_vector(stationary_states[i],
//                              std::string("grond_state_") +
//                                  Utilities::int_to_string(i));
//   data_out.build_patches();
//
//   fs::path file_path = results_path / "ground_states.vtu";
//   std::ofstream output(file_path);
//   data_out.write_vtu(output);
// }


/*------------------ Explicit templates instantiation -------------------*/

template class Model<1>;
template class Model<2>;
template class Model<3>;

// vim: ts=2 sts=2 sw=2 fdc=2
