#ifndef TIME_DEPENDENT_HEADER
#define TIME_DEPENDENT_HEADER

#include "model.hpp"
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparse_matrix.h>
#include <filesystem>

using std::complex;
using namespace dealii;

template <int dim>
class TDSE : public BaseProblem
{
public:
  TDSE(Model<dim> &model);
  void run() override;
  void setup_system();
  void assemble_matrices();
  void do_time_step();
  void output_results() const;
  void save_potential() const;

private:
  Model<dim> &model;

  const dealii::ParameterHandler &parameters;
  std::filesystem::path results_path;

  dealii::Triangulation<dim> &mesh;
  const dealii::FiniteElement<dim> &fe;
  dealii::DoFHandler<dim> &dof_handler;

  const std::vector<dealii::Vector<double>> &initial_states;

  dealii::AffineConstraints<complex<double>> constraints;

  dealii::SparsityPattern sparsity_pattern;

  dealii::SparseMatrix<complex<double>> matrix_next;
  dealii::SparseMatrix<complex<double>> matrix_current;

  dealii::Vector<complex<double>> solution;
  dealii::Vector<complex<double>> system_rhs;

  double       time;
  const double time_step;
  unsigned int timestep_number;
  const double end_time;
};

#endif // TIME_DEPENDENT_HEADER
// vim: ts=2 sts=2 sw=2
