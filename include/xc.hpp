#ifndef XC_HEADER
#define XC_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <deal.II/lac/vector.h>

namespace xc {

  /**
   * @brief Zero-frequency limit of longitudinal fxc.
   *
   * Returns @p n*f0, where
   * @verbatim f0 = d^2(n*exc(n))/dn^2 @endverbatim
   * given by Eq.(i) in
   * Ref. E. K. U. Gross and W. Kohn, Phys. Rev. Lett. 55, 2850 (1985).
   *
   * @param n13 Qubit root of electron density, @p n^{1/3}
   */
  template <class Number>
  dealii::Vector<Number> get_nf0 (dealii::Vector<Number> n13);

}

#endif // XC_HEADER
// vim: ts=2 sts=2 sw=2
