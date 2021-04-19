#include "ground_state.hpp"
#include "model.hpp"

#include <deal.II/base/function_parser.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/slepc_solver.h>

// #include <deal.II/base/logstream.h>
// #include <deal.II/base/quadrature_lib.h>

// #include <deal.II/base/function.h>

// #include <deal.II/base/utilities.h>
// #include <deal.II/grid/tria.h>
// #include <deal.II/grid/grid_in.h>
// #include <deal.II/grid/grid_generator.h>
// #include <deal.II/grid/tria_accessor.h>
// #include <deal.II/grid/tria_iterator.h>
// #include <deal.II/dofs/dof_handler.h>
// #include <deal.II/dofs/dof_accessor.h>

// #include <deal.II/fe/fe_q.h>
// #include <deal.II/fe/fe_values.h>
// #include <deal.II/numerics/matrix_tools.h>
// #include <deal.II/lac/affine_constraints.h>
// #include <deal.II/lac/full_matrix.h>
// #include <deal.II/base/index_set.h>
// #include <deal.II/lac/petsc_sparse_matrix.h>
// #include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <fstream>
#include <filesystem>

using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;

template <int dim>
EigenvalueProblem<dim>::EigenvalueProblem (Model<dim> &model)
  : model(model),
    parameters(model.parameters),
    mesh(model.mesh),
    fe(model.get_fe()),
    dof_handler(model.dof_handler)
{}

template <int dim>
void EigenvalueProblem<dim>::setup_system() //{{{
{
  cout << "Global mesh refinement steps: "
       << parameters.get_integer("Global mesh refinement steps") << endl;
  mesh.refine_global( parameters.get_integer("Global mesh refinement steps") );
  dof_handler.distribute_dofs(fe);

  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  stiffness_matrix.reinit(dof_handler.n_dofs(),
                          dof_handler.n_dofs(),
                          dof_handler.max_couplings_between_dofs());
  mass_matrix.reinit(dof_handler.n_dofs(),
                     dof_handler.n_dofs(),
                     dof_handler.max_couplings_between_dofs());

  IndexSet eigenfunction_index_set = dof_handler.locally_owned_dofs();
  eigenfunctions.resize(
    parameters.get_integer("Number of eigenvalues/eigenfunctions"));
  for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
    eigenfunctions[i].reinit( eigenfunction_index_set, MPI_COMM_WORLD );

  eigenvalues.resize( eigenfunctions.size() );

} //}}}

template <int dim>
void EigenvalueProblem<dim>::assemble_system() //{{{
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  FunctionParser<dim> potential;
  potential.initialize(FunctionParser<dim>::default_variable_names(),
                       parameters.get("Potential"),
                       typename FunctionParser<dim>::ConstMap());

  std::vector<double> potential_values(n_q_points);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_stiffness_matrix = 0;
      cell_mass_matrix      = 0;

      potential.value_list(fe_values.get_quadrature_points(),
                           potential_values);

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            cell_stiffness_matrix(i, j) +=
                (
                  1. / 2 *                            // 1/2
                  fe_values.shape_grad(i, q_point) *  // grad phi_i
                  fe_values.shape_grad(j, q_point)    // grad phi_j
                  +
                  potential_values[q_point] *         // V(r)
                  fe_values.shape_value(i, q_point) * // phi_i
                  fe_values.shape_value(j, q_point)   // phi_j
                )
                * fe_values.JxW(q_point);             // J*dr

            cell_mass_matrix(i, j) +=
                (
                  fe_values.shape_value(i, q_point) * // grad(phi_i)
                  fe_values.shape_value(j, q_point)   // grad(phi_j)
                )
                * fe_values.JxW(q_point);             // J*dr
          }

      // Now that we have the local matrix contributions, we transfer them
      // into the global objects and take care of zero boundary constraints:
      cell->get_dof_indices(local_dof_indices);

      constraints.distribute_local_to_global(cell_stiffness_matrix,
                                             local_dof_indices,
                                             stiffness_matrix);
      constraints.distribute_local_to_global(cell_mass_matrix,
                                             local_dof_indices,
                                             mass_matrix);
    }

  // At the end of the function, we tell PETSc that the matrices have now
  // been fully assembled and that the sparse matrix representation can now
  // be compressed as no more entries will be added:
  stiffness_matrix.compress(VectorOperation::add);
  mass_matrix.compress(VectorOperation::add);


  double min_spurious_eigenvalue = std::numeric_limits<double>::max(),
         max_spurious_eigenvalue = -std::numeric_limits<double>::max();

  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    if (constraints.is_constrained(i))
      {
        const double ev         = stiffness_matrix(i, i) / mass_matrix(i, i);
        min_spurious_eigenvalue = std::min(min_spurious_eigenvalue, ev);
        max_spurious_eigenvalue = std::max(max_spurious_eigenvalue, ev);
      }

  std::cout << "   Spurious eigenvalues are all in the interval "
            << "[" << min_spurious_eigenvalue << ","
            << max_spurious_eigenvalue << "]" << std::endl;
} //}}}

template <int dim>
unsigned int EigenvalueProblem<dim>::solve() //{{{
{
  // Assigning convergence control
  SolverControl                    solver_control(dof_handler.n_dofs(), 1e-9);
  SLEPcWrappers::SolverKrylovSchur eigensolver(solver_control);

  eigensolver.set_which_eigenpairs(EPS_SMALLEST_REAL);

  eigensolver.set_problem_type(EPS_GHEP);

  eigensolver.solve(stiffness_matrix,
                    mass_matrix,
                    eigenvalues,
                    eigenfunctions,
                    eigenfunctions.size());

  for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
    eigenfunctions[i] /= eigenfunctions[i].linfty_norm();

  // Return the number of iterations it took to converge
  return solver_control.last_step();
} //}}}

template <int dim>
void EigenvalueProblem<dim>::output_results() const //{{{
{
  // Save eigenfunctions and interpolated potential.
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);

  for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
    data_out.add_data_vector(eigenfunctions[i],
                             std::string("eigenfunction_") +
                                 Utilities::int_to_string(i));

  Vector<double> projected_potential(dof_handler.n_dofs());
  {
    FunctionParser<dim> potential;
    potential.initialize(FunctionParser<dim>::default_variable_names(),
                         parameters.get("Potential"),
                         typename FunctionParser<dim>::ConstMap());
    VectorTools::interpolate(dof_handler, potential, projected_potential);
  }
  data_out.add_data_vector(projected_potential, "interpolated_potential");

  data_out.build_patches();

  // Write to file
  fs::path file_path = model.results_path / "eigenvectors.vtu";
  std::ofstream output(file_path);
  data_out.write_vtu(output);
} //}}}

template <int dim>
void EigenvalueProblem<dim>::run() //{{{
{
  setup_system();

  cout << "   Number of active cells:       "
       << mesh.n_active_cells() << endl
       << "   Number of degrees of freedom: " << dof_handler.n_dofs()
       << endl;

  assemble_system();

  const unsigned int n_iterations = solve();
  cout << "   Solver converged in " << n_iterations << " iterations."
       << endl;

  // Convert files from dealii::PETScWrappers::MPI::Vector to dealii::Vector.
  {
    model.ground_states.resize( eigenfunctions.size() );
    for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
      model.ground_states[i] = eigenfunctions[i];

    // model.ground_states.clear();
    // model.ground_states.reserve( eigenfunctions.size() );
    // for (auto eigfun : eigenfunctions)
    //   model.ground_states.emplace_back( eigfun );
  }

  output_results();

  cout << endl;
  for (unsigned int i = 0; i < eigenvalues.size(); ++i)
    cout << "      Eigenvalue " << i << " : " << eigenvalues[i] << endl;
  cout << endl;

} //}}}

// template <int dim>
// void EigenvalueProblem<dim>::test()
// {}

// Explicit template instantiation
template class EigenvalueProblem<1>;
template class EigenvalueProblem<2>;
template class EigenvalueProblem<3>;


// vim: ts=2 sts=2 sw=2
