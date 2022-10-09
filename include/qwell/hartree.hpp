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
class Hartree : public HartreeBase
{
// private:
  // virtual void setup_boundary_constraints() override;
  virtual void setup_boundary_conditions() override;
};

#endif // HARTREE_HEADER
// vim: ts=2 sts=2 sw=2
