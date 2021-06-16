#ifndef KOHN_SHAM_HEADER
#define KOHN_SHAM_HEADER

#include "models/base_model.hpp"

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table_handler.h>
#include <filesystem>
// #include <utility>


struct KohnShamOrbitals
{
  std::vector<double> eigenvalues;
  std::vector <dealii::Vector<double>> wavefunctions;
  std::array <double, 2> spurious_eigenvalues_interval;
};


template <int dim>
class KohnSham : public BaseProblem
{
public:
  struct Parameters
  {
    unsigned int number_of_electrons;
  };

  KohnSham (Model<dim> &model,
            const Parameters & parameters,
            const dealii::Vector<double> & fe_potential,
            const dealii::Function<dim> & fun_potential);

  KohnSham (Model<dim> & model,
            const Parameters & parameters,
            const dealii::Vector<double> & fe_potential,
            const dealii::Function<dim> * fun_potential);

  KohnSham (Model<dim> &model,
            const Parameters & parameters,
            const dealii::Vector<double> & fe_potential);

  KohnShamOrbitals run();

// private:
  void setup_system();
  void assemble_system();
  unsigned int solve();
  std::vector<dealii::Vector<double>> get_wavefunctions() const;

  Model<dim> & model;
  const Parameters & parameters;

  /// External potential in the form of the finite element function.
  const dealii::Vector<double> * fe_potential;

  /// External potential in the form of the @p dealii::Function<dim> object.
  const dealii::Function<dim> * fun_potential;

  dealii::Triangulation<dim> & mesh;
  const dealii::FiniteElement<dim> & fe;
  dealii::DoFHandler<dim> & dof_handler;

  dealii::AffineConstraints<double> constraints;

  dealii::PETScWrappers::SparseMatrix stiffness_matrix, mass_matrix;
  std::vector <dealii::PETScWrappers::MPI::Vector> eigenfunctions;
  std::vector<double> eigenvalues;

  double min_spurious_eigenvalue;
  double max_spurious_eigenvalue;

  dealii::TableHandler output_table;
};


#endif // KOHN_SHAM_HEADER
// vim: ts=2 sts=2 sw=2
