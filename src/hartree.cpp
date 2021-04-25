#include "hartree.hpp"
#include "model.hpp"

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <math.h>
#include <iostream>
#include <fstream>
#include <filesystem>

using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;

template <int dim>
HartreeProblem<dim>::HartreeProblem (Model<dim> &model)
  : model(model),
    parameters(model.parameters),
    mesh(model.mesh),
    fe(model.get_fe()),
    dof_handler(model.dof_handler),
    stationary_states(model.stationary_states),
    number_of_electrons(stationary_states.size())
{}



template <int dim>
void HartreeProblem<dim>::setup_system()
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
void HartreeProblem<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients | update_JxW_values);

  MatrixCreator::create_laplace_matrix(dof_handler,
                                       quadrature_formula,
                                       system_matrix,
                                       (const Function<dim> * const)nullptr,
                                       constraints);

  // const unsigned int dofs_per_cell = fe_values.dofs_per_cell; // fe.dofs_per_cell;
  //
  //FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
  //Vector<double>     cell_rhs (dofs_per_cell);
  //
  // std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  //
  // for (const auto &cell : dof_handler.active_cell_iterators())
  // {
  //   fe_values.reinit(cell);
  //   cell_matrix = 0;
  //   cell_rhs    = 0;
  //
  //   // Assemble cell matrix
  //   for (const unsigned int q_index : fe_values.quadrature_point_indices())
  //     for (const unsigned int i : fe_values.dof_indices())
  //       for (const unsigned int j : fe_values.dof_indices())
  //       {
  //         cell_matrix(i, j) +=
  //             fe_values.shape_grad(i, q_index) * // grad phi_i(r_q)
  //             fe_values.shape_grad(j, q_index) * // grad phi_j(r_q)
  //             fe_values.JxW(q_index);            // J*dr
  //       }
  //
  //   cell->get_dof_indices(local_dof_indices);
  //   constraints.distribute_local_to_global(cell_matrix,
  //                                          local_dof_indices,
  //                                          system_matrix);
  // }; // end loop over active cells


  // Setup the right hand side of the system.
  system_rhs.reinit( stationary_states[0].size() );

  for (unsigned int i = 0; i < system_rhs.size(); ++i)
    for (unsigned int state = 0; state < number_of_electrons; ++state)
    {
      // 4 * pi * density
      system_rhs[i] += 4 * M_PI * stationary_states[state][i] *
                                  stationary_states[state][i];
    }
  constraints.condense(system_rhs);

}



template <int dim>
void HartreeProblem<dim>::solve()
{
  SolverControl solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}



template <int dim>
void HartreeProblem<dim>::output_results() const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "Hartree_potential");
  data_out.build_patches();

  // Write to file
  fs::path file_path = model.results_path / "Hartree_potential.vtu";
  std::ofstream output(file_path);
  data_out.write_vtu(output);
}



template <int dim>
void HartreeProblem<dim>::run()
{
  setup_system();
  assemble_system();
  solve();
  output_results();
}



// template <int dim>
// void HartreeProblem<dim>::test()
// {}



// Explicit template instantiation
template class HartreeProblem<1>;
template class HartreeProblem<2>;
template class HartreeProblem<3>;


// vim: ts=2 sts=2 sw=2
