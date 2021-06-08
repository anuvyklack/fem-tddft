#include "hartree.hpp"
#include "model.hpp"

#include <deal.II/fe/fe_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
// #include <deal.II/lac/solver_cg.h>
// #include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_direct.h>

#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <filesystem>

using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;


template <int dim>
Hartree<dim>::Hartree (Model<dim> & model, dealii::Vector<double> & density)
  : model(model),
    mesh(model.mesh),
    model_fe(model.get_fe()),
    dof_handler(mesh),
    density(density)
{
  // Make local finite element object be the same type as in 'model',
  // but of the 2 time higher polinomial order.
  fe = FETools::get_fe_by_name<dim>(
           model.parameters.fe_type +
           "(" +  std::to_string(model.parameters.fe_order * 2) + ")"
       );

  dof_handler.distribute_dofs(*fe);
}



template <int dim>
void Hartree<dim>::setup_system()
{
  DynamicSparsityPattern dsp( dof_handler.n_dofs() );
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  // Dirichlet boundary conditions.
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();
}



template <int dim>
void Hartree<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(fe->degree + 1);

  FEValues<dim> fe_values(
      *fe, quadrature_formula,
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
  FETools::interpolate(model.dof_handler, density,
                       dof_handler, system_rhs);
  system_rhs *= 4 * M_PI;

  constraints.condense(system_rhs);
}



template <int dim>
void Hartree<dim>::solve()
{
  // SolverControl solver_control(1000, 1e-12);
  // SolverCG<Vector<double>> solver(solver_control);
  // solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  // Direct solver:
  // --------------
  SparseDirectUMFPACK direct_solver;
  direct_solver.solve(system_matrix, system_rhs);
  solution = system_rhs;

  constraints.distribute(solution);
}



template <int dim>
dealii::Vector<double>
Hartree<dim>::run()
{
  setup_system();
  assemble_system();
  solve();

  // We don't need 'system_rhs' vector any more, so to not allocate memory for
  // a temporary vector,  we use 'system_rhs' vector as temporary storage during
  // interpolation 'solution' vector back to 'model.dof_handler':

  // resize the vector
  // system_rhs.reinit(0); // free memory but take time
  system_rhs.reinit(model.dof_handler.n_dofs());

  FETools::interpolate(dof_handler, solution,
                       model.dof_handler, system_rhs);
  return system_rhs;
}



/*------------------ Explicit templates instantiation -------------------*/

template class Hartree<1>;
template class Hartree<2>;
template class Hartree<3>;

// vim: ts=2 sts=2 sw=2
