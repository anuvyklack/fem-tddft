#ifndef EXTERNAL_POTENTIALS_HEADER
#define EXTERNAL_POTENTIALS_HEADER

#include <deal.II/base/function.h>


template <int dim>
class PointCharge : public dealii::Function<dim>
{
public:
  PointCharge (double charge_magnitude,
               dealii::Point<dim> charge_location);

  virtual double value (const dealii::Point<dim> & p,
                        const unsigned int component=0) const override;

private:
  double charge_magnitude;
  dealii::Point<dim> charge_location;
};


#endif // EXTERNAL_POTENTIALS_HEADER
// vim: ts=2 sts=2 sw=2
