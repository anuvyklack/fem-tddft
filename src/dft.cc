#include "dft.hpp"
#include "hartree.hpp"
#include "utilities.hpp"

#include <deal.II/lac/solver_control.h>
// #include <deal.II/numerics/time_dependent.h>
#include <iostream>

using namespace dealii;
using std::cin, std::cout, std::endl;

DFT_Parameters::DFT_Parameters()
  : ParameterAcceptor("/DFT")
{
  add_parameter(
      "Number of electrons",
      number_of_electrons,
      "How many electrons (i.e. the number of eigenvalues/eigenfunctions) "
      "to be computed.");

  parse_parameters_call_back.connect( [&](){initialized = true;} );
}


template <int dim>
DFT<dim>::DFT (Model<dim> & model, DFT_Parameters & prm)
  : model(model), parameters(prm)
{
  AssertThrow(prm.initialized, ExcNotInitialized())

  // Set hartree potential and density vector sizes.
  // And set there initial values to 0.
  hartree_potential.reinit( model.dof_handler.n_dofs() );
  density.reinit( model.dof_handler.n_dofs());
}



template <int dim>
DFT<dim>::DFT (Model<dim> & model, DFT_Parameters & prm,
               const dealii::Function<dim> & external_potential)
  : model(model), parameters(prm), external_potential(&external_potential)
{
  AssertThrow(prm.initialized, ExcNotInitialized())

  // Set hartree potential and density vectors sizes,
  // and set there initial values to 0.
  hartree_potential.reinit( model.dof_handler.n_dofs() );
  density.reinit( model.dof_handler.n_dofs());
}



template <int dim>
void DFT<dim>::run()
{
  SolverControl solver_control(100, 1e-9);

  dealii::Vector<double> delta;
  delta.reinit(hartree_potential);

  unsigned int it = 0; // the number of iteration
  SolverControl::State solver_state = SolverControl::iterate;
  while (solver_state == SolverControl::iterate)
    {
      ++it;
      delta = hartree_potential;

      cout
      << "   Number of active cells:       " << model.mesh.n_active_cells() << endl
      << "   Number of degrees of freedom: " << model.dof_handler.n_dofs() << endl;

      solve_Kohn_Sham_problem();
      solve_Hartree_problem();

      delta -= hartree_potential;
      // double residual = delta.linfty_norm();
      double residual = delta.l2_norm();

      solver_state = solver_control.check(it, residual);

      if (solver_state == SolverControl::failure)
        AssertThrow(false, SolverControl::NoConvergence(it, residual));
    };

  // // Output the interval all spurious eigenvalues are in.
  // cout << "   Spurious eigenvalues are all in the interval "
  //      << "[" << kohn_sham_output.spurious_eigenvalues_interval.first << ","
  //             << kohn_sham_output.spurious_eigenvalues_interval.second << "]"
  //      << endl;
}



template <int dim>
void DFT<dim>::solve_Kohn_Sham_problem()
{
  typename KohnSham<dim>::Parameters ksprm {parameters.number_of_electrons};
  KohnSham<dim> ks (model, ksprm, hartree_potential, external_potential);

  kohn_sham_orbitals = ks.run();
  // data.density = get_density();
}



template <int dim>
void DFT<dim>::solve_Hartree_problem()
{
  // Mix old and new densities.
  temp.reinit(density);
  temp = density;
  density = get_density();
  temp = mixer(density, temp, 0.5);

  Hartree<dim> hartree {model, parameters, temp};
  hartree_potential = hartree.run();
}



template <int dim>
dealii::Vector<double>
DFT<dim>::get_density() const
{
  // Number of degrees of freedom, i.e. vector size.
  unsigned int n_dofs = model.dof_handler.n_dofs();
  dealii::Vector<double> density (n_dofs);

  for (unsigned int i = 0; i < n_dofs; ++i)
    for (unsigned int wf = 0; wf < parameters.number_of_electrons; ++wf)
      density[i] += kohn_sham_orbitals.wavefunctions[wf][i] *
                    kohn_sham_orbitals.wavefunctions[wf][i];

  return density;
}


/*--------------- Explicit templates instantiation --------------------------*/

template class DFT<1>;
template class DFT<2>;
template class DFT<3>;

// vim: ts=2 sts=2 sw=2
