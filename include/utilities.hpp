#ifndef UTILITIES_HEADER
#define UTILITIES_HEADER

#include <deal.II/fe/mapping_q1_eulerian.h>
#include <magic_enum.hpp>

/// @brief Check if this app is running in parellel with "mpirun".
/// Exit the application with error if true.
void exit_if_mpi (char** &argv);


/*--------------- Tlmplate and inline functions -----------------------------*/

/**
 * @brief Pulay Mixing
 *
 * Returns
 * @code
 *   theta * new_vec + (1-theta) * old_vec
 * @endcode
 */
template <typename Number>
dealii::Vector<Number>
mixer (const dealii::Vector<Number> & new_vec,
       const dealii::Vector<Number> & old_vec,
       double theta)
{
  dealii::Vector<Number> result (new_vec.size());
  result.add(   theta, new_vec );
  result.add( 1-theta, old_vec );
  return result;
}



/**
 * @brief Return a string suitable to pass to @p dealii::Patterns::Selection()
 * constructor.
 *
 * From enumeration
 * @code
 *   enum MyEnum {first, second, third};
 * @endcode
 * construct the string of the form
 * @code
 *   "first|second|third"
 * @endcode
 * which is suitable to pass to @link dealii::Patterns::Selection @endlink
 * constructor.
 *
 * @tparam ENUM the enumeration to process.
 */
template <typename Enum>
std::string enum_to_selection_string()
{
  std::string result;
  for (auto name : magic_enum::enum_names<Enum>())
    result += static_cast<std::string>(name) + '|';
  if (result.back() == '|') result.pop_back(); // Delete the last '|" symbol.
  return result;
}


#endif // UTILITIES_HEADER
// vim: ts=2 sts=2 sw=2
