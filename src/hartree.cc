#include "hartree.hpp"
#include "model.hpp"

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <filesystem>

using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;


template <int dim>
Hartree<dim>::Hartree (Model<dim> & model,
                       const DFT_Parameters & parameters,
                       DFT_Data & data)
  : model(model),
    parameters(parameters),
    mesh(model.mesh),
    fe(model.get_fe()),
    dof_handler(model.dof_handler),
    density(data.density)
{}



template <int dim>
void Hartree<dim>::setup_system()
{
  DynamicSparsityPattern dsp( dof_handler.n_dofs() );
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  // Dirichlet boundary conditions
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();
}



template <int dim>
void Hartree<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(
      fe, quadrature_formula,
      update_values | update_gradients | update_JxW_values
  );

  MatrixCreator::create_laplace_matrix(
      dof_handler,
      quadrature_formula,
      system_matrix,
      (const Function<dim> * const)nullptr,
      constraints
  );

  // Setup the right hand side of the system.
  system_rhs.reinit(density);
  system_rhs = density;
  system_rhs *= 4 * M_PI;

  constraints.condense(system_rhs);
}



template <int dim>
void Hartree<dim>::solve()
{
  SolverControl solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  constraints.distribute(solution);
}



template <int dim>
dealii::Vector<double>
Hartree<dim>::run()
{
  setup_system();
  assemble_system();
  solve();
  return solution;
}



/*--------------- Explicit templates instantiation --------------------------*/

template class Hartree<1>;
template class Hartree<2>;
template class Hartree<3>;

// vim: ts=2 sts=2 sw=2
