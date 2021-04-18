#ifndef TIME_DEPENDENT_HPP
#define TIME_DEPENDENT_HPP

#include "model.hpp"
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparse_matrix.h>

using std::complex;
using namespace dealii;

template <int dim>
class TDSE : BaseProblem
{
public:
  TDSE(Model<dim> &model);
  void run() override;
  void restore_data();  // Read data from archive.
  void setup_system();
  void assemble_matrices();
  void do_time_step();
  void output_results() const;
  void save_potential() const;

private:
  const std::string results_folder = "results";

  ParameterHandler parameters;

  Triangulation<dim> mesh;
  FE_Q<dim> fe;
  // FiniteElement<dim> *fe;
  DoFHandler<dim> dof_handler {mesh};
  std::vector<dealii::Vector<double>> initial_states;

  AffineConstraints<complex<double>> constraints;

  SparsityPattern sparsity_pattern;

  SparseMatrix<complex<double>> matrix_next;
  SparseMatrix<complex<double>> matrix_current;

  dealii::Vector<complex<double>> solution;
  dealii::Vector<complex<double>> system_rhs;

  double time = 0;
  const double time_step {0.00001};
  unsigned int timestep_number {0};
  const double end_time;
};

#endif // TIME_DEPENDENT_HPP

// vim: ts=2 sts=2 sw=2
