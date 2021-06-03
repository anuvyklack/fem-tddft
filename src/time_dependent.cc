#include "time_dependent.hpp"

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

#include <iostream>
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;
using std::cin, std::cout, std::endl;
using std::complex;
using namespace dealii;



template <int dim>
class Potential : public Function<dim>
{
  virtual double value(const Point<dim> & p,
                       const unsigned int component=0) const override
  {
    // Does nothing, but suppress the compiler 
    // warning about the unused variable.
    (void)component;
    (void)p;

    Assert(component == 0, ExcIndexRange(component, 0, 1));

    // The distance between origin and point.
    // return 1. / Point<dim>().distance(p);
    return 0;
  };
};



template <int dim>
TDSE<dim>::TDSE(Model<dim> & model)
  : model(model),
    parameters(model.parameters),
    mesh(model.mesh),
    fe(model.get_fe()),
    dof_handler(model.dof_handler),
    initial_states(model.stationary_states),
    time(0),
    time_step(1./64),
    timestep_number(0),
    end_time(1)
{
  results_path = model.results_path / "time_dependent_solution",
  fs::create_directory(results_path);
}



template <int dim>
void TDSE<dim>::run()
{
  setup_system();
  assemble_matrices();

  for (unsigned int i = 0; i < solution.size(); ++i){
    solution[i] = initial_states[0][i];
  };
  // VectorTools::interpolate(dof_handler, InitialValues<dim>(), solution);
  output_results();

  // double time = 0;
  for (; time <= end_time; time += time_step)
  {
    ++timestep_number;
    cout << "Time step " << timestep_number << " at t=" << time << endl;
    do_time_step();

    // if (timestep_number % 1 == 0)
      output_results();
  }
}



template <int dim>
void TDSE<dim>::setup_system()
{
  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  matrix_next.reinit(sparsity_pattern);
  matrix_current.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  // Dirichlet boundary conditions
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();
}



template <int dim>
void TDSE<dim>::assemble_matrices()
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell; // fe.dofs_per_cell;

  FullMatrix<complex<double>> cell_matrix_lhs {dofs_per_cell, dofs_per_cell};
  FullMatrix<complex<double>> cell_matrix_rhs {dofs_per_cell, dofs_per_cell};

  std::vector<double> potential_values (fe_values.n_quadrature_points);
  const Potential<dim> potential;

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    cell_matrix_lhs = complex<double>(0.);
    cell_matrix_rhs = complex<double>(0.);

    fe_values.reinit(cell);

    potential.value_list(fe_values.get_quadrature_points(),
                         potential_values);

    const std::complex<double> i {0, 1}; // imaginary unit

    for (const unsigned int q_index : fe_values.quadrature_point_indices())
      for (const unsigned int k : fe_values.dof_indices())
        for (const unsigned int l : fe_values.dof_indices())
        {

          cell_matrix_lhs(k, l) +=
              (
                fe_values.shape_value(k, q_index) *   // phi_k
                fe_values.shape_value(l, q_index)     // phi_l
                +
                i * time_step *                       // i*dt
                (
                  1./2 *                              // 1/2
                  fe_values.shape_grad(k, q_index) *  // grad phi_k
                  fe_values.shape_grad(l, q_index)    // grad phi_l
                  +
                  potential_values[q_index] *         // V(r)
                  fe_values.shape_value(k, q_index) * // phi_k
                  fe_values.shape_value(l, q_index)   // phi_l
                )
              ) * fe_values.JxW(q_index);             // J*dr


          cell_matrix_rhs(k, l) +=
              (
                fe_values.shape_value(k, q_index) *   // phi_k
                fe_values.shape_value(l, q_index)     // phi_l
                -
                i * time_step *                       // i*dt
                (
                  1./2 *                              // 1/2
                  fe_values.shape_grad(k, q_index) *  // grad phi_k
                  fe_values.shape_grad(l, q_index)    // grad phi_l
                  +
                  potential_values[q_index] *         // V(r)
                  fe_values.shape_value(k, q_index) * // phi_k
                  fe_values.shape_value(l, q_index)   // phi_l
                )
              ) * fe_values.JxW(q_index);             // J*dr

        }

    cell->get_dof_indices(local_dof_indices);
    constraints.distribute_local_to_global(cell_matrix_lhs,
                                           local_dof_indices,
                                           matrix_next);
    constraints.distribute_local_to_global(cell_matrix_rhs,
                                           local_dof_indices,
                                           matrix_current);
  } // end loop over active cells
}



template <int dim>
void TDSE<dim>::do_time_step()
{
                    // OutVector,  InVector
  matrix_current.vmult(system_rhs, solution);

  // SparseDirectUMFPACK takes the matrix and the right hand side vector and
  // returns the solution in the same vector it takes.
  SparseDirectUMFPACK direct_solver;
  direct_solver.solve(matrix_next, system_rhs);

  solution = system_rhs;
}



template <int dim>
void TDSE<dim>::save_potential() const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  Vector<double> projected_potential(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, Potential<dim>{}, projected_potential);
  data_out.add_data_vector(projected_potential, "interpolated_potential");
  data_out.build_patches();

  fs::path file_path = fs::current_path() / results_path / "potential.vtk";
  std::ofstream output(file_path);
  data_out.write_vtk(output);
}



namespace DataPostprocessors
{
  template <int dim>
  class ProbabilityDensity : public DataPostprocessorScalar<dim>
  {
  public:
    ProbabilityDensity()
        : DataPostprocessorScalar<dim>("Density", update_values)
    {};

    virtual void evaluate_vector_field (
        const DataPostprocessorInputs::Vector<dim> &inputs,
        std::vector<Vector<double>> &computed_quantities
    ) const override
    {
      Assert(computed_quantities.size() == inputs.solution_values.size(),
             ExcDimensionMismatch(computed_quantities.size(),
                                  inputs.solution_values.size()));

      for (unsigned int q = 0; q < computed_quantities.size(); ++q)
      {
        Assert(computed_quantities[q].size() == 1,
               ExcDimensionMismatch(computed_quantities[q].size(), 1));
        Assert(inputs.solution_values[q].size() == 2,
               ExcDimensionMismatch(inputs.solution_values[q].size(), 2));

        const std::complex<double> psi(inputs.solution_values[q](0),
                                       inputs.solution_values[q](1));
        computed_quantities[q](0) = std::norm(psi);
      }
    };
  }; // class ProbabilityDensity
} // namespace DataPostprocessors



template <int dim>
void TDSE<dim>::output_results() const
{
  const DataPostprocessors::ProbabilityDensity<dim> density;
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "Psi");
  data_out.add_data_vector(solution, density);
  data_out.build_patches();
  data_out.set_flags(DataOutBase::VtkFlags(time, timestep_number));

  const std::string filename =
      "td_solution_" + Utilities::int_to_string(timestep_number, 3) + ".vtu";
  fs::path file_path =
      fs::current_path() / results_path / filename;
  std::ofstream output(file_path);

  data_out.write_vtu(output);
}



// Explicit template instantiation
template class TDSE<1>;
template class TDSE<2>;
template class TDSE<3>;

// vim: ts=2 sts=2 sw=2
