#ifndef QUANTUM_WELL_HEADER
#define QUANTUM_WELL_HEADER

#include "models/base_model.hpp"
#include "atomic_units.hpp"
#include "parameters_parsing.hpp"
#include "double_output_stream.hpp"

#include <deal.II/base/function_parser.h>


template <int dim>
class QuantumWell : public Model<dim>
{
public:
  /// The Meyers Singleton.
  class Parameters : public dealii::ParameterAcceptor
  {
  public:
    static Parameters & get_parameters()
    {
      static QuantumWell<dim>::Parameters prm;
      return prm;
    }

    double width = 400;    ///< The width of the well.

    double eps = 13.0;   ///< Dielectric constant.
    double meff = 0.067; ///< Effective mass in units m0.

    double sheet_density = 1.0e11;  ///< Sheet electron density.

    std::string seed_density_fun = "2 * (ns/width) * cos(pi * x / width)^2";

    dealii::FunctionParser<dim> seed_density;
    std::map<std::string, double> seed_density_constants;


    std::unique_ptr<EffectiveAU> effau;

    bool initialized = false;

  private:
    Parameters(); /* Constructor should be private! */
    Parameters (const Parameters& other) = delete;
    Parameters& operator=(const Parameters&) = delete;

    virtual void parse_parameters (dealii::ParameterHandler & prm) override;
  };

  QuantumWell (typename Model<dim>::Parameters & model_parameters,
               Parameters & parameters);

  // ~QuantumWell();

  const Parameters & parameters;
  EffectiveAU & effau;  ///< Effective atomic units handler.

  // using Model<dim>::results_path;
  using Model<dim>::out;
  using Model<dim>::mesh;
  using Model<dim>::dof_handler;
  using Model<dim>::fe_ptr;

private:
  virtual void set_mesh () override;
};


#endif // QUANTUM_WELL_HEADER
// vim: ts=2 sts=2 sw=2
