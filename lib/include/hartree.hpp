#ifndef HARTREE_HEADER
#define HARTREE_HEADER

#include "models/base_model.hpp"
#include "dft.hpp"

#include <deal.II/grid/tria.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparse_matrix.h>

#include <filesystem>

template <int dim>
class Hartree : public BaseProblem
{
public:
  Hartree (Model<dim> & model, dealii::Vector<double> & density);

  dealii::Vector<double> run();

// private:
  void setup_system();
  void assemble_system();
  void solve();

  Model<dim> & model;

  dealii::Triangulation<dim> & mesh;
  const dealii::FiniteElement<dim> & model_fe;

  /**
   * The local for this class finite element object should be the same type
   * as in @link Model @endlink class object given into constructor, but
   * of the 2 time higher order.
   */
  std::unique_ptr <dealii::FiniteElement<dim>> fe;

  dealii::DoFHandler<dim> dof_handler;

  dealii::Vector<double> & density;

  dealii::AffineConstraints<double> constraints;

  dealii::SparsityPattern      sparsity_pattern;
  dealii::SparseMatrix<double> system_matrix;

  dealii::Vector<double> solution;
  dealii::Vector<double> system_rhs;  ///< 4 * pi * density.

};

#endif // HARTREE_HEADER
// vim: ts=2 sts=2 sw=2
