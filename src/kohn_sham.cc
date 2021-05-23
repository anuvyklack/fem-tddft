#include "kohn_sham.hpp"
#include "model.hpp"
#include "dft.hpp"

#include <deal.II/base/function_parser.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/slepc_solver.h>

#include <iostream>
#include <fstream>
#include <filesystem>

using namespace dealii;
using std::cin, std::cout, std::endl;
namespace fs = std::filesystem;



template <int dim>
KohnSham<dim>::KohnSham (Model<dim> & model,
                         const DFT_Parameters & parameters,
                         DFT_Data & data)

  : model(model),
    parameters(parameters),
    hartree_potential(data.hartree_potential),
    mesh(model.mesh),
    fe(model.get_fe()),
    dof_handler(model.dof_handler)
{}



template <int dim>
void KohnSham<dim>::setup_system()
{
  // Dirichlet boundary conditions.
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  stiffness_matrix.reinit(dof_handler.n_dofs(),
                          dof_handler.n_dofs(),
                          dof_handler.max_couplings_between_dofs());

  mass_matrix.reinit(dof_handler.n_dofs(),
                     dof_handler.n_dofs(),
                     dof_handler.max_couplings_between_dofs());

  IndexSet eigenfunction_index_set = dof_handler.locally_owned_dofs();
  eigenfunctions.resize( parameters.number_of_electrons );
  for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
    eigenfunctions[i].reinit( eigenfunction_index_set, MPI_COMM_WORLD );

  eigenvalues.resize( eigenfunctions.size() );
}



template <int dim>
void KohnSham<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell; // fe.dofs_per_cell;

  FullMatrix<double> cell_stiffness_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_mass_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  // FunctionParser<dim> potential;
  // potential.initialize(FunctionParser<dim>::default_variable_names(),
  //                      parameters.get("Potential"),
  //                      typename FunctionParser<dim>::ConstMap());
  //
  // std::vector<double> potential_values (fe_values.n_quadrature_points);

  std::vector<double> hartree_values (fe_values.n_quadrature_points);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    cell_stiffness_matrix = 0;
    cell_mass_matrix      = 0;

    // potential.value_list(fe_values.get_quadrature_points(),
    //                      potential_values);

    // Get only currently needed local part of hartree patential.
    fe_values.get_function_values( hartree_potential, hartree_values );

    for (const unsigned int q_index : fe_values.quadrature_point_indices())
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
        {
          cell_stiffness_matrix(i, j) +=
              (
                1./2 *                              // 1/2
                fe_values.shape_grad(i, q_index) *  // grad phi_i
                fe_values.shape_grad(j, q_index)    // grad phi_j
                +
                hartree_values[q_index] *           // V_Hartree(r)
                fe_values.shape_value(i, q_index) * // phi_i
                fe_values.shape_value(j, q_index)   // phi_j
                // +
                // potential_values[q_index] *         // V(r)
                // fe_values.shape_value(i, q_index) * // phi_i
                // fe_values.shape_value(j, q_index)   // phi_j
              )
              * fe_values.JxW(q_index);             // J*dr

          cell_mass_matrix(i, j) +=
              (
                fe_values.shape_value(i, q_index) * // grad(phi_i)
                fe_values.shape_value(j, q_index)   // grad(phi_j)
              )
              * fe_values.JxW(q_index);             // J*dr
        }

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
        const double ev = stiffness_matrix(i,i) / mass_matrix(i,i);
        min_spurious_eigenvalue = std::min(min_spurious_eigenvalue, ev);
        max_spurious_eigenvalue = std::max(max_spurious_eigenvalue, ev);
      }

  cout
  << "   Spurious eigenvalues are all in the interval "
  << "[" << min_spurious_eigenvalue << "," << max_spurious_eigenvalue << "]"
  << std::endl;
}



template <int dim>
unsigned int KohnSham<dim>::solve()
{
  // Assigning convergence control.
  SolverControl solver_control(dof_handler.n_dofs(), 1e-9);
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
}



template <int dim>
std::vector<dealii::Vector<double>>
KohnSham<dim>::get_orbitals() const
{
  std::vector<dealii::Vector<double>> orbitals;

  // Convert eigenfunctions from
  // 'dealii::PETScWrappers::MPI::Vector'
  // to 'dealii::Vector'.
  {
    orbitals.resize( eigenfunctions.size() );
    for (unsigned int i=0; i < eigenfunctions.size(); ++i)
      orbitals[i] = eigenfunctions[i];

    // Second version:

    // orbitals.reserve( eigenfunctions.size() );
    // for (auto eigfun : eigenfunctions)
    //   orbitals.emplace_back( eigfun );
  }

  // Normalize result wave functions.
  for (dealii::Vector<double> & state : orbitals)
    state /= state.l2_norm();

  return orbitals;
}



template <int dim>
std::vector<dealii::Vector<double>>
KohnSham<dim>::run()
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

  cout << endl;
  for (unsigned int i = 0; i < eigenvalues.size(); ++i)
    cout << "      Eigenvalue " << i << " : " << eigenvalues[i] << endl;
  cout << endl;

  return get_orbitals();
}



/*--------------- Explicit templates instantiation --------------------------*/

template class KohnSham<1>;
template class KohnSham<2>;
template class KohnSham<3>;


// vim: ts=2 sts=2 sw=2
