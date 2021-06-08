#ifndef DFT_HEADER
#define DFT_HEADER

#include "model.hpp"
#include "double_output_stream.hpp"
#include "kohn_sham.hpp"
// #include "hartree.hpp"
#include "parameters_parsing.hpp"
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/convergence_table.h>

/**
 * @brief The parameters for DFT calculation.
 *
 * The Meyers Singleton.
 */
class DFT_Parameters : public dealii::ParameterAcceptor
{
public:
  static DFT_Parameters & get_parameters()
  {
    static DFT_Parameters prm;
    return prm;
  }

  unsigned int number_of_electrons = 2;

  bool initialized = false;

private:
  DFT_Parameters();
  DFT_Parameters (const DFT_Parameters& other) = delete;
  DFT_Parameters& operator= (const DFT_Parameters&) = delete;
};



// template <int dim>
// struct DFT_Data
// {
//   std::vector<dealii::Vector<double>> kohn_sham_orbitals;
//   dealii::Vector<double> density;
//   dealii::Vector<double> hartree_potential;
//   const dealii::Function<dim>* external_potential {nullptr};
// };



template <int dim>
class DFT : public BaseProblem
{
public:
  DFT (Model<dim> & model, DFT_Parameters & parameters);

  DFT (Model<dim> & model, DFT_Parameters & prm,
       const dealii::Function<dim> & external_potential);

  void run();

  KohnShamOrbitals kohn_sham_orbitals;
  dealii::Vector<double> density;
  dealii::Vector<double> hartree_potential;

// private:
  void solve_Kohn_Sham_problem();
  void solve_Hartree_problem();
  void refine_mesh();

  dealii::Vector<double> get_density() const;

  /// Return full potential energy.
  dealii::Vector<double> get_potential() const;

  Model<dim> & model;
  const DFT_Parameters & parameters;
  const dealii::Function<dim>* external_potential {nullptr};

  double_ostream & out;

  /**
   * Auxiliary vector for some additioanl calculations.
   * We create it here to avoid creation and consequently memory allocation
   * (which is time-consuming) on each iteration when we need it.
   */
  dealii::Vector<double>  temp;

  dealii::ConvergenceTable convergence_table;
};



// template <int dim>
// class DFT_Adaptive : public DFT
// {
//
// };


#endif // DFT_HEADER
// vim: ts=2 sts=2 sw=2
