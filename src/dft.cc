#include "dft.hpp"
#include "kohn_sham.hpp"
#include "hartree.hpp"
#include "output_results.hpp"

#include <deal.II/lac/solver_control.h>
// #include <deal.II/numerics/time_dependent.h>

using namespace dealii;

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

  // Set initial hartree potential to 0.
  data.hartree_potential.reinit( model.dof_handler.n_dofs() );
  data.hartree_potential = 0;
}



template <int dim>
void DFT<dim>::run()
{
  SolverControl solver_control(100, 1e-9);

  dealii::Vector<double> delta;
  delta.reinit(data.hartree_potential);

  unsigned int it = 0; // the number of iteration
  SolverControl::State solver_state = SolverControl::iterate;
  while (solver_state == SolverControl::iterate)
    {
      ++it;
      delta = data.hartree_potential;

      solve_Kohn_Sham_problem();
      calculate_density();
      solve_Hartree_problem();

      delta -= data.hartree_potential;
      double residual = delta.linfty_norm();

      solver_state = solver_control.check(it, residual);

      if (solver_state == SolverControl::failure)
        AssertThrow(false, SolverControl::NoConvergence(it, residual));
    };
}



template <int dim>
void DFT<dim>::solve_Kohn_Sham_problem()
{
  KohnSham<dim> ks {model, parameters, data};
  data.kohn_sham_orbitals = ks.run();
}



template <int dim>
void DFT<dim>::solve_Hartree_problem()
{
  Hartree<dim> hartree {model, parameters, data};
  data.hartree_potential =  hartree.run();
}



template <int dim>
void DFT<dim>::calculate_density()
{
  unsigned int n_dofs = model.dof_handler.n_dofs();
  data.density.reinit(n_dofs);

  for (unsigned int i = 0; i < n_dofs; ++i)
    for (unsigned int e_num = 0;
         e_num < parameters.number_of_electrons; ++e_num)
      data.density[i] += data.kohn_sham_orbitals[e_num][i] *
                         data.kohn_sham_orbitals[e_num][i];
}


/*--------------- Explicit templates instantiation --------------------------*/

template class DFT<1>;
template class DFT<2>;
template class DFT<3>;

// vim: ts=2 sts=2 sw=2
