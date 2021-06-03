#include "xc.hpp"
#include <math.h>

using namespace dealii;

const double ax = -(3./4.) * cbrt(3. / M_PI);
const double ac = (log(2.) - 1.) / (2. * M_PI * M_PI);
const double B1 = 21.7392245 * cbrt(4. * M_PI / 3.);
const double B2 = 20.4562557 *  cbrt(16. * M_PI * M_PI / 9.);


template <class Number>
Vector<Number> get_nf0 (Vector<Number> n13)
{
  Vector<Number> result (n13.size());
  for(unsigned int m = 0; m < n13.size(); ++m)
    {
      double cr1 = n13(m);
      double cr2 = cr1 * cr1;
      double cr3 = cr1 * cr2;

      // exchange
      result(m) = 4. * ax * cr1 / 9.;

      // correlation
      double arg = 1. + B1*cr1 + B2*cr2;
      result(m) += ac * cr1
                   * (4. * B1
                      + (3. * B1 * B1 + 10.*B2) * cr1
                      + 10. * B1 * B2 * cr2
                      + 6. * B2 * B2 * cr3)
                   / (9. * arg * arg);
    }
  return result;
}

/*--------------- Explicit templates instantiation ----------------------*/

template Vector<double> get_nf0 (Vector<double>);

// vim: ts=2 sts=2 sw=2
