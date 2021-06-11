#ifndef DFT_HEADER
#define DFT_HEADER

#include "models/base_model.hpp"
#include "double_output_stream.hpp"
#include "kohn_sham.hpp"
// #include "hartree.hpp"
#include "parameters_parsing.hpp"
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/convergence_table.h>

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
  /**
   * @brief The parameters for DFT calculation.
   */
  class Parameters : public dealii::ParameterAcceptor
  {
  public:
    Parameters();
    Parameters (const Parameters& other) = delete;
    Parameters& operator= (const Parameters&) = delete;
  
    unsigned int number_of_electrons = 2;
    unsigned int max_convergence_steps = 100;
  
    bool initialized = false;
  };


  DFT (Model<dim> & model, Parameters & parameters);

  DFT (Model<dim> & model, Parameters & prm,
       const dealii::Function<dim> & external_potential);

  DFT (Model<dim> & model, Parameters & prm,
       const dealii::Function<dim> *external_potential,
       const dealii::Function<dim> *seed_density);

  void run();

  KohnShamOrbitals kohn_sham_orbitals;
  dealii::Vector<double> density;
  dealii::Vector<double> hartree_potential;

// private:
  void solve_Kohn_Sham_problem();
  void solve_Hartree_problem();
  void refine_mesh();

  dealii::Vector<double> get_density() const;

  dealii::Vector<double> get_hartree_plus_xc_potential();

  Model<dim> & model;
  const Parameters & parameters;
  const dealii::Function<dim> * external_potential {nullptr};
  bool use_seed_density = false;

  double_ostream & out;

  dealii::ConvergenceTable convergence_table;

  /**
   * Auxiliary vector for some additioanl calculations.
   * We create it here to avoid creation and consequently memory allocation
   * (which is time-consuming) on each iteration when we need it.
   */
  dealii::Vector<double> temp;
};



// template <int dim>
// class DFT_Adaptive : public DFT
// {
//
// };


#endif // DFT_HEADER
// vim: ts=2 sts=2 sw=2
