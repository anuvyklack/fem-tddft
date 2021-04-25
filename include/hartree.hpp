#ifndef HARTREE_HPP
#define HARTREE_HPP

#include "model.hpp"

#include <deal.II/grid/tria.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparse_matrix.h>

#include <filesystem>

template <int dim>
class HartreeProblem : public BaseProblem
{
public:
  HartreeProblem (Model<dim> &model);
  void run() override;

  // void test();

private:
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;

  Model<dim> &model;
  const dealii::ParameterHandler &parameters;

  dealii::Triangulation<dim> &mesh;
  const dealii::FiniteElement<dim> &fe;
  dealii::DoFHandler<dim> &dof_handler;

  const std::vector<dealii::Vector<double>> &stationary_states;

  dealii::AffineConstraints<double> constraints;

  const unsigned int number_of_electrons;

  dealii::SparsityPattern      sparsity_pattern;
  dealii::SparseMatrix<double> system_matrix;

  dealii::Vector<double> solution;
  dealii::Vector<double> system_rhs;  ///< Actually, this is 4 * pi * density.

};

#endif
// vim: ts=2 sts=2 sw=2
