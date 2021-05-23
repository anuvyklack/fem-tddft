#ifndef KOHN_SHAM_HEADER
#define KOHN_SHAM_HEADER

#include "model.hpp"
#include "dft.hpp"

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/base/parameter_handler.h>
#include <filesystem>

template <int dim>
class KohnSham : public BaseProblem
{
public:
  KohnSham (Model<dim> &model,
            const DFT_Parameters & parameters,
            DFT_Data & data);

  std::vector<dealii::Vector<double>> run();

// private:
  void setup_system();
  void assemble_system();
  unsigned int solve();
  std::vector<dealii::Vector<double>> get_orbitals() const;

  Model<dim> & model;
  const DFT_Parameters & parameters;
  dealii::Vector<double> & hartree_potential;

  dealii::Triangulation<dim> & mesh;
  const dealii::FiniteElement<dim> & fe;
  dealii::DoFHandler<dim> & dof_handler;

  dealii::AffineConstraints<double> constraints;

  dealii::PETScWrappers::SparseMatrix stiffness_matrix, mass_matrix;
  std::vector <dealii::PETScWrappers::MPI::Vector> eigenfunctions;
  std::vector<double> eigenvalues;
};


#endif // KOHN_SHAM_HEADER
// vim: ts=2 sts=2 sw=2
