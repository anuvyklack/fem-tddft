#ifndef XC_HEADER
#define XC_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <deal.II/lac/vector.h>

namespace xc {

  /**
   * @brief LDA exchange-correlation potential.
   *
   * @note The cubic root of the electron density @p cr=pow(rho,1/3.)
   * is passed as an argument due to performance considerations.
   *
   * @param n13 Qubit root of electron density, @p n^{1/3}
   */
  template <class Number>
  dealii::Vector<Number> get_VxcLDA (const dealii::Vector<Number> n13);

} // namespace xc

#endif // XC_HEADER
// vim: ts=2 sts=2 sw=2
