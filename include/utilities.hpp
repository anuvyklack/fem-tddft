#ifndef UTILITIES_HEADER
#define UTILITIES_HEADER

#include "model.hpp"
#include <magic_enum.hpp>

/// @brief Check if this app is running in parellel with "mpirun".
/// Exit the application with error if true.
void exit_if_mpi (char** &argv);

/*--------------- Tlmplate and inline functions -----------------------------*/

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
template <typename ENUM>
std::string enum_to_selection_string()
{
  std::string result;
  for (auto name : magic_enum::enum_names<ENUM>())
    result += static_cast<std::string>(name) + '|';
  if (result.back() == '|') result.pop_back(); // Delete the last '|" symbol.
  return result;
}


#endif // UTILITIES_HEADER
// vim: ts=2 sts=2 sw=2
