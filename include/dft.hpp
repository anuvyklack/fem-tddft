#ifndef DFT_HEADER
#define DFT_HEADER

#include "model.hpp"
#include <deal.II/base/parameter_acceptor.h>
#include "parameters_parsing.hpp"

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



struct DFT_Data
{
  std::vector<dealii::Vector<double>> kohn_sham_orbitals;
  dealii::Vector<double> density;
  dealii::Vector<double> hartree_potential;
};



template <int dim>
class DFT : public BaseProblem
{
public:
  DFT (Model<dim> & model, DFT_Parameters & parameters);
  void run();

// private:
  void solve_Kohn_Sham_problem();
  void solve_Hartree_problem();
  void calculate_density();
  void refine_mesh();

  Model<dim> & model;
  const DFT_Parameters & parameters;
  DFT_Data data;
};



// template <int dim>
// class DFT_Adaptive : public DFT
// {
//
// };


#endif // DFT_HEADER
// vim: ts=2 sts=2 sw=2
