#include "dft.hpp"
#include "hartree.hpp"
#include "xc.hpp"
#include "utilities.hpp"

#include <deal.II/lac/solver_control.h>
// #include <deal.II/numerics/time_dependent.h>
#include <iostream>
#include <cmath>

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
  add_parameter(
      "Maximum DFT convergence steps",
      max_convergence_steps,
      "The maximum number of iterations the solver is allowed to do trying "
      "to achieve the solution convergence.");

  parse_parameters_call_back.connect( [&](){initialized = true;} );
}


template <int dim>
DFT<dim>::DFT (Model<dim> & model, DFT_Parameters & prm)
  : model(model), parameters(prm), out(model.out)
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
  : model(model),
    parameters(prm),
    external_potential(&external_potential),
    out(model.out)
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
  out << "=============================================" << std::endl;

  SolverControl solver_control(parameters.max_convergence_steps, 1e-9);

  dealii::Vector<double> delta;
  delta.reinit(hartree_potential);

  // For convenience, we define next variables:
  double n_cells = model.mesh.n_active_cells();
  double n_dofs = model.dof_handler.n_dofs();
  double L2_error, Linfty_error;

  unsigned int iteration = 0; // the number of iteration
  SolverControl::State solver_state = SolverControl::iterate;
  while (solver_state == SolverControl::iterate)
    {
      ++iteration;

      out << endl
          << "Iteration " << iteration << ':' << endl
          << "   Number of active cells:       " << n_cells << endl
          << "   Number of degrees of freedom: " << n_dofs  << endl;

      delta = hartree_potential;

      out << "   Solving Kohn-Sham equation." << endl;
      solve_Kohn_Sham_problem();

      out << "   Solving Poisson's equation." << endl;
      solve_Hartree_problem();

      delta -= hartree_potential;
      Linfty_error = delta.linfty_norm();
      L2_error = delta.l2_norm();

      solver_state = solver_control.check(iteration, L2_error);

      // if (solver_state == SolverControl::failure)
      //   AssertThrow(false, SolverControl::NoConvergence(iteration, residual));

      convergence_table.add_value("#", iteration);
      convergence_table.add_value("cells", n_cells);
      convergence_table.add_value("dofs", n_dofs);
      convergence_table.add_value("L2", L2_error);
      convergence_table.add_value("Linfty", Linfty_error);
    };

  out << endl << "=============================================" << endl;

  // convergence_table.add_column_to_supercolumn("L2", "error");
  // convergence_table.add_column_to_supercolumn("Linfty", "error");

  convergence_table.set_precision("cells", 0);
  convergence_table.set_precision("dofs", 0);
  convergence_table.set_precision("L2", 3);
  convergence_table.set_precision("Linfty", 3);

  convergence_table.set_scientific("L2", true);
  convergence_table.set_scientific("Linfty", true);

  convergence_table.set_tex_caption("cells", "\\# cells");
  convergence_table.set_tex_caption("dofs", "\\# dofs");
  convergence_table.set_tex_caption("L2", "@f$L^2@f$-error");
  convergence_table.set_tex_caption("Linfty", "@f$L^\\infty@f$-error");

  convergence_table.set_tex_format("cells", "r");
  convergence_table.set_tex_format("dofs", "r");

  convergence_table.evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate);
  convergence_table.evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate_log2);

  // convergence_table.evaluate_convergence_rates("Linfty", ConvergenceTable::reduction_rate);
  // convergence_table.evaluate_convergence_rates("Linfty", ConvergenceTable::reduction_rate_log2);

  out << endl;
  convergence_table.write_text(out.out1, TableHandler::org_mode_table);
  convergence_table.write_text(out.out2, TableHandler::org_mode_table);


  out << endl << "=============================================" << endl << endl;

  for (unsigned int i = 0; i < kohn_sham_orbitals.eigenvalues.size(); ++i)
    out << "Kohn-Sham eigenvalue " << i << " : "
         << kohn_sham_orbitals.eigenvalues[i] << endl;
  out << endl;

  out << "Spurious eigenvalues are all in the interval "
       << "["
       << kohn_sham_orbitals.spurious_eigenvalues_interval[0] << ", "
       << kohn_sham_orbitals.spurious_eigenvalues_interval[1]
       << "]"
       << endl;
}



template <int dim>
void DFT<dim>::solve_Kohn_Sham_problem()
{
  temp = get_hartree_plus_xc_potential();
  typename KohnSham<dim>::Parameters ksprm {parameters.number_of_electrons};
  KohnSham<dim> ks (model, ksprm, temp, external_potential);

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

  Hartree<dim> hartree {model, temp};
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



template <int dim>
dealii::Vector<double>
DFT<dim>::get_hartree_plus_xc_potential()
{
  unsigned int n_dofs = model.dof_handler.n_dofs();

  // Calculate the qubic root of the density.
  temp.reinit(density);
  for (unsigned int i = 0; i < n_dofs; ++i)
    temp[i] = cbrt(density[i]);

  return xc::get_VxcLDA(temp);
}



/*------------------ Explicit templates instantiation -------------------*/

template class DFT<1>;
template class DFT<2>;
template class DFT<3>;

// vim: ts=2 sts=2 sw=2
