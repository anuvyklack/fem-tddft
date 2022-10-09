#ifndef ATOMIC_UNITS_HEADER
#define ATOMIC_UNITS_HEADER
// -------------------------------------------------------
// Copyright (C) 2019 by Dmitry R. Gulevich
// You may use or modify this code but not distribute it.
// -------------------------------------------------------
#include <iostream>
#include <cmath>

// Values of physical constants from NIST (https://physics.nist.gov)
const double BOHR_in_Angstrom = 0.52917721067; ///< Bohr radius in Angstrom
const double HARTREE_in_eV = 27.21138602; ///< eV
const double AUOT_in_sec = 2.418884326509e-17; ///< atomic unit of time in seconds
// const double eps0 = 8.854187817e-12; // electric constant, F/m

/**
 * @brief Effective atomic units
 */
struct EffectiveAU
{
  EffectiveAU (double eps, double meff)
    : eps(eps), m_eff(meff)
  {}

  const double eps;   ///< Dielectric constant,
  const double m_eff; ///< Effective mass,

  /// Effective charge in units of e.
  const double e_eff = 1. / std::sqrt(eps);

  /// Effective Bohr radius in Angstrom.
  const double a0 = BOHR_in_Angstrom / (m_eff * e_eff * e_eff);

  /// Effective Bohr radius in 'cm'.
  const double a0_cm = a0 * 1.e-8;

  /// Effective Hartree energy in 'meV'.
  const double Eh = 1.e3 * HARTREE_in_eV * m_eff * std::pow(e_eff,4);

  /// V/m, effective atomic unit of electric field by Eh/e0a0 (note e not eeff!).
  const double E_au = 1.e-3 * Eh / (a0 * 1.e-10);

  const double E_mV_nm = 1.e-6 * E_au; /// mV/nm

  /// Unit time in fs.
  const double ut_fs = 1.e15 * AUOT_in_sec * HARTREE_in_eV / (Eh * 1.e-3);


  /// Convenient conversion from physical units to effective a.u.
  double to_au (double value, std::string units)
  {
    if      (units == "meV")      return value / Eh;
    else if (units == "Angstrom") return value / a0;
    else if (units == "cm")       return value / a0_cm;
    else if (units == "mV/nm")    return value / E_mV_nm;
    else
    {
      std::cout << "Error. Supplied units are not implemented." << std::endl;
      throw std::invalid_argument( "Supplied units are not implemented." );
      return 0.;
    }
  }


  /// Convenient conversion from from effective a.u. to physical units.
  double from_au (double value, std::string units)
  {
    if      (units == "meV")      return value * Eh;
    else if (units == "Angstrom") return value * a0;
    else if (units == "cm")       return value * a0_cm;
    else if (units == "mV/nm")    return value * E_mV_nm;
    else
    {
      std::cout << "Error. Supplied units are not implemented." << std::endl;
      throw std::invalid_argument( "Supplied units are not implemented." );
      return 0.;
    }
  }

}; // struct EffectiveAU

#endif // ATOMIC_UNITS_HEADER
// vim: ts=2 sts=2 sw=2
