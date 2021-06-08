#include "xc.hpp"
#include <math.h>

using namespace dealii;

const double ax = -(3./4.) * cbrt(3. / M_PI);
const double ac = (log(2.) - 1.) / (2. * M_PI * M_PI);
const double B1 = 21.7392245 * cbrt(4. * M_PI / 3.);
const double B2 = 20.4562557 *  cbrt(16. * M_PI * M_PI / 9.);


// LDA xc potential
template <class Number>
Vector<Number> get_VxcLDA (Vector<Number> n13) 
{
  Vector<Number> result (n13.size());
  for (unsigned int m = 0; m < n13.size(); ++m) 
    {
      // exchange
      double cr1 = n13[m]; // cubic root
      result[m] = 4 * ax * cr1 / 3.;

      // correlation
      double cr2 = cr1*cr1; // cubic root squared
      double arg = 1. + B1*cr1 + B2*cr2;
      result[m] += ac * ((B1*cr1 + 2*B2*cr2) / (3.*arg) + log(arg));
    }
  return result;
}


/*------------------ Explicit templates instantiation -------------------*/

template Vector<double> get_VxcLDA (Vector<double>);

// vim: ts=2 sts=2 sw=2
