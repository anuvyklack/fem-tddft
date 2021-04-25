#ifndef STATIONARY_SCHRODINGER_HPP
#define STATIONARY_SCHRODINGER_HPP

#include "model.hpp"

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/base/parameter_handler.h>
#include <filesystem>

template <int dim>
class EigenvalueProblem : public BaseProblem
{
public:
  EigenvalueProblem (Model<dim> &model);
  void run() override;

  // void test();

private:
  void setup_system();
  void assemble_system();
  unsigned int solve();
  void output_results() const;

  Model<dim> &model;

  const dealii::ParameterHandler &parameters;

  dealii::Triangulation<dim> &mesh;
  const dealii::FiniteElement<dim> &fe;
  dealii::DoFHandler<dim> &dof_handler;

  dealii::AffineConstraints<double> constraints;

  dealii::PETScWrappers::SparseMatrix stiffness_matrix, mass_matrix;
  std::vector <dealii::PETScWrappers::MPI::Vector> eigenfunctions;
  std::vector<double> eigenvalues;
};

#endif // STATIONARY_SCHRODINGER_HPP

// vim: ts=2 sts=2 sw=2
