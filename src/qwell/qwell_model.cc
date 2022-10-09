#include "models/quantum_well.hpp"
#include <deal.II/grid/grid_generator.h>

using namespace dealii;
using std::endl;


template <int dim>
QuantumWell<dim>::Parameters::Parameters()
  : ParameterAcceptor("/Quantum Well")
{
  add_parameter(
      "Quantum well width",
      width,
      "The width of the quantum well in angstroms.");
  add_parameter(
      "Dielectric constant",
      eps);
  add_parameter(
      "Effective mass",
      meff,
      "Effective mass in units m0.");
  add_parameter(
      "Sheet electron density",
      sheet_density);
  add_parameter(
       "Seed density",
       seed_density,
       "A functional description of the initial guess of density "
       "for DFT calculation.");

  parse_parameters_call_back.connect( [&](){initialized = true;} );
}



template <int dim>
void QuantumWell<dim>::Parameters::parse_parameters (ParameterHandler & prm)
{
  (void)prm;

  // Setup effective atomic units.
  effau = std::make_unique<EffectiveAU>(eps, meff);

  width = effau->to_au(width, "Angstrom");
  sheet_density = sheet_density * effau->a0_cm * effau->a0_cm;
}



template<int dim>
QuantumWell<dim>::QuantumWell (
    typename Model<dim>::Parameters & model_parameters,
    Parameters & parameters
)
  : Model<dim>(model_parameters),
    parameters(parameters),
    effau(*parameters.effau)
{
  AssertThrow(parameters.initialized, ExcNotInitialized())

  set_mesh();
  out << "Global mesh refinement steps: "
      << std::to_string( Model<dim>::parameters.global_mesh_refinement_steps )
      << endl << endl;
  mesh.refine_global( Model<dim>::parameters.global_mesh_refinement_steps );

  dof_handler.distribute_dofs(*fe_ptr);

  // Parse seed density function.
  {
    std::map<std::string, double> seed_density_constants;

    seed_density_constants["pi"] = numbers::PI;
    seed_density_constants["ns"] = parameters.sheet_density;
    seed_density_constants["width"] = parameters.width;

    seed_density.initialize(FunctionParser<dim>::default_variable_names(),
                            parameters.seed_density,
                            seed_density_constants);
  }
}



template <int dim>
void QuantumWell<dim>::set_mesh ()
{
  GridGenerator::hyper_cube(mesh, 0, parameters.width,
                            true /* colorize boundaries */);
}



/*--------------------- Explicit templates instantiation ---------------------*/

template class QuantumWell<1>;
template class QuantumWell<2>;
template class QuantumWell<3>;

// vim: ts=2 sts=2 sw=2 fdc=2
