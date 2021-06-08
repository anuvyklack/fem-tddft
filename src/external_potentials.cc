#include "external_potentials.hpp"

using namespace dealii;



template <int dim>
PointCharge<dim>::PointCharge (double charge_magnitude,
                               dealii::Point<dim> charge_location)
  : charge_magnitude(charge_magnitude),
    charge_location(charge_location)
{}



template <int dim>
double PointCharge<dim>::value(const Point<dim> & p,
                               const unsigned int component) const
{
  (void)component; /* Does nothing, but suppress the compiler
                      warning about the unused variable. */
  Assert(component == 0, ExcIndexRange(component, 0, 1));

  return charge_magnitude / charge_location.distance(p);
}


/*------------------ Explicit templates instantiation -------------------*/

template class PointCharge<2>;
template class PointCharge<3>;
