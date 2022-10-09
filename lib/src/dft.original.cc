#include "dft.hpp"
#include "hartree.hpp"
#include "xc.hpp"
#include "utilities.hpp"

#include <deal.II/lac/solver_control.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/sparse_direct.h>
// #include <deal.II/numerics/time_dependent.h>
#include <iostream>
#include <cmath>

using namespace dealii;
using std::cin, std::cout, std::endl;


template <int dim>
DFT<dim>::Parameters::Parameters()
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

  // add_parameter(
  //     "Mixer",
  //     mixer_type,
  //     "The type of density mixing to use. "
  //     "Possible values:\n"
  //     "\t- pulay\n"
  //     "\t- linear");

  enter_subsection("Pulay density mixing");
    add_parameter(
        "Mixer size",
        density_max_size,
        "The number of last densities to use in Pulay mixing to calculate "
        "the new density.");
    add_parameter(
        "Linear mixing steps",
        linear_mixing_steps,
        "The number of iterations at the beginning of calculation during which "
        "the linear mixing will be used. Could not be less than 1.");
  leave_subsection();


  parse_parameters_call_back.connect( [&](){initialized = true;} );
}



// template <int dim>
// void DFT<dim>::Parameters::parse_parameters (ParameterHandler & prm)
// {
//   (void)prm;
//   if (mixer_type)
// }



template <int dim>
DFT<dim>::DFT (Model<dim> & model, Parameters & prm)
  : DFT(model, prm, nullptr, nullptr)
{}



template <int dim>
DFT<dim>::DFT (Model<dim> & model, Parameters & prm,
               const dealii::Function<dim> &external_potential)
  : DFT(model, prm, &external_potential, nullptr)
{}



template <int dim>
DFT<dim>::DFT (Model<dim> & model, Parameters & prm,
               const Function<dim> *external_potential,
               const Function<dim> *seed_density)
  : model(model),
    parameters(prm),
    external_potential(external_potential),
    out(model.out)
{
  AssertThrow(prm.initialized, ExcNotInitialized())

  // Set hartree potential and density vectors size,
  // and set their initial values to 0.
  hartree_potential.reinit( model.dof_handler.n_dofs() );

  temp.reinit( model.dof_handler.n_dofs() );
  if (seed_density)
    {
      use_seed_density = true;
      VectorTools::interpolate(model.dof_handler, *seed_density, temp);
      density.push_back(temp);
    }
}



template <int dim>
void DFT<dim>::run()
{
  out << "=============================================" << endl;

  SolverControl solver_control(parameters.max_convergence_steps, 1e-9);

  dealii::Vector<double> delta;
  delta.reinit(hartree_potential);

  // For convenience, we define next variables:
  double n_cells = model.mesh.n_active_cells();
  double n_dofs = model.dof_handler.n_dofs();
  double L2_error, Linfty_error;

  unsigned int iteration = 0; // the number of iteration

  if (use_seed_density)
    {
      out << endl
          << "Iteration " << iteration << ':' << endl
          << "   Use initial density." << endl
          << "   Solving Poisson's equation." << endl;

      Hartree<dim> hartree {model, density.back()};
      hartree_potential = hartree.run();
    }

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

      calculate_density();

      out << "   Solving Poisson's equation." << endl;
      solve_Hartree_problem();

      delta -= hartree_potential;
      Linfty_error = delta.linfty_norm();
      L2_error = delta.l2_norm();

      out << "   Error:                        " << L2_error << endl;

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
  // convergence_table.evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate_log2);

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
}



template <int dim>
void DFT<dim>::solve_Hartree_problem()
{
  temp = pulay_mixer();
  // temp = linear_mixer();
  Hartree<dim> hartree {model, temp};
  hartree_potential = hartree.run();
}



template <int dim>
const dealii::Vector<double>&
DFT<dim>::get_density() const
{
  return density.back();
}



template <int dim>
void DFT<dim>::calculate_density()
{
  // Number of degrees of freedom, i.e. vector size.
  unsigned int n_dofs = model.dof_handler.n_dofs();
  temp.reinit(n_dofs);

  for (unsigned int i = 0; i < n_dofs; ++i)
    for (unsigned int wf = 0; wf < parameters.number_of_electrons; ++wf)
      temp[i] += kohn_sham_orbitals.wavefunctions[wf][i] *
                 kohn_sham_orbitals.wavefunctions[wf][i];

  density.push_back(temp);

  // If the length of the density array is greater than 1, we calculate
  // the difference between current density and previous.
  if (density.size() > 1)
    {
      temp -= density[density.size()-2];
      density_error.push_back(temp);
    }

  // Keep density array size no more than density_max_size.
  if (density.size() > density_max_size)
    {
      density.pop_front();
      Assert(density.size() == density_max_size,
             ExcDimensionMismatch(density.size(), density_max_size));
    }

  // Make density and density_error arrays be the equal length.
  if (density_error.size() > density.size())
    density_error.pop_front();
  else if (density.size() > density_error.size())
    density.pop_front();

  // Check if we obtained desired.
  Assert(density.size() == density_error.size(),
         ExcDimensionMismatch(density.size(), density_error.size()));

  update_mixer_matrix();
}



template <int dim>
void DFT<dim>::update_mixer_matrix()
{
  if (density.size() == 1)
    return;
  else if (density.size() == 2)
    {
      // Ensure that density and density_error arrays are of the same length.
      Assert(density.size() == density_error.size(),
             ExcDimensionMismatch(density.size(), density_error.size()));

      mixer_matrix.reinit(3,3);

      unsigned int max_index = mixer_matrix.m() - 1;

      // Set diagonal elements.
      for (unsigned int i=0; i<2; ++i)
          mixer_matrix(i,i) = density_error[i] * density_error[i];

      // Set non-diagonal elements.
      for (unsigned int i = 0; i < 2; ++i)
        for (unsigned int j = i+1; j < 2; ++j)
          {
            mixer_matrix(i,j) = density_error[i] * density_error[j];
            mixer_matrix(j,i) = mixer_matrix(i,j);
          }
      for (unsigned int i = 0; i < max_index; ++i)
        {
          mixer_matrix(max_index, i) = 1;
          mixer_matrix(i, max_index) = 1;
        }

      // cout << endl;
      // mixer_matrix.print_formatted(cout, 1, true, 8);
      // cout << endl;

      return;
    }
  // else if (density.size() < density_max_size)
  else if (density.size() == mixer_matrix.m())
    {
      // Ensure that density and density_error arrays are of the same length.
      Assert(density.size() == density_error.size(),
             ExcDimensionMismatch(density.size(), density_error.size()));

      // Ensure that matrix is square.
      Assert(mixer_matrix.m() == mixer_matrix.n(),
             ExcDimensionMismatch(mixer_matrix.m(), mixer_matrix.n()));

      // Assert(density.size() == mixer_matrix.m(),
      //        ExcDimensionMismatch(density.size(), mixer_matrix.m()));

      mixer_matrix.grow_or_shrink( density.size() + 1 );

      unsigned int max_index = mixer_matrix.m() - 1;

      for (unsigned int i=0; i < max_index; ++i)
        {
          mixer_matrix(max_index, i) = 1;
          mixer_matrix(i, max_index) = 1;
        }
    }
  else
    {
      Assert(density.size() == density_max_size,
             ExcDimensionMismatch(density.size(), density_max_size));

      // mixer_matrix.fill(mixer_matrix,
      //                   0, 0,  // offsets of the input matrix
      //                   1, 1); // offsets of the being copied matrix

      unsigned int max_index = mixer_matrix.m() - 1;
      for (unsigned int i=1; i < max_index; ++i)
        for (unsigned int j=1; j < max_index; ++j)
          mixer_matrix(i-1,j-1) = mixer_matrix(i,j);
    }


  unsigned int max_index = mixer_matrix.m() - 1;

  mixer_matrix(max_index-1, max_index-1) =
      density_error.back() * density_error.back();

  for (unsigned int i=0; i < max_index; ++i)
    {
      mixer_matrix(max_index-1, i) = density_error[i] * density_error.back();
      mixer_matrix(i, max_index-1) = mixer_matrix(max_index-1, i);
    }

   // cout << endl;
   // mixer_matrix.print_formatted(cout, 1, true, 8);
   // cout << endl;
}



template <int dim>
dealii::Vector<double>
DFT<dim>::pulay_mixer()
{
  // if (density.size() == 1) return density.back();
  if (density.size() < linear_mixing_steps)
    return linear_mixer();

  cout << "   Use Pulay mixing." << endl;

  // Right hand side vector.
  dealii::Vector<double> c (mixer_matrix.m()); // coefficients
  c[c.size()-1] = 1;

  LAPACKFullMatrix< double > lu_matrix = mixer_matrix;
  lu_matrix.compute_lu_factorization();
  lu_matrix.solve(c);

  c.grow_or_shrink( c.size() - 1 );

  // cout << "Coeffitients: " << c.size() << endl;
  // cout << "Density: " << density.size() << endl;
  // cout << "Density error: " << density_error.size() << endl;

  temp.reinit( model.dof_handler.n_dofs() );
  for (unsigned int i = 0; i < density.size(); ++i)
    temp.add(c[i], density[i]);

  // c.print(cout);
  // density[0].print(cout);
  // density[1].print(cout);
  // temp.print(cout);

  return temp;
}



template <int dim>
dealii::Vector<double>
DFT<dim>::linear_mixer()
{
  cout << "   Use linear mixing." << endl;

  if (density.size() == 1)
    return density.back();
  else
    {
      temp.reinit( model.dof_handler.n_dofs() );
      double theta = 0.5;
      temp.add(     theta, density.back() );
      temp.add( 1 - theta, density[density.size()-2] );
      return temp;
    }
}



template <int dim>
dealii::Vector<double>
DFT<dim>::get_hartree_plus_xc_potential()
{
  unsigned int n_dofs = model.dof_handler.n_dofs();

  // Calculate the qubic root of the density.
  temp.reinit(density.back());
  for (unsigned int i = 0; i < n_dofs; ++i)
    temp[i] = cbrt(density.back()[i]);

  temp = xc::get_VxcLDA(temp);
  temp += hartree_potential;
  return temp;
}

/*------------------ Explicit templates instantiation -------------------*/

template class DFT<1>;
template class DFT<2>;
template class DFT<3>;

// vim: ts=2 sts=2 sw=2
