# ${CMAKE_BINARY_DIR}
# ${PROJECT_SOURCE_DIR}

# PROJECT(schrodinger VERSION 0.0.3)
PROJECT(fem-tddft)

SET(LIBRARY femdft)
SET(TARGETS qwell test)

SET(DEAL_II_DIR ~/software/dealii-9.3.0/build)


# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 3.12.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Global output directories for targets.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Search for deal.II library.
# ---------------------------
FIND_PACKAGE(deal.II 9.3.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

# Is deal.II correctly configured?
# --------------------------------
IF(NOT DEAL_II_WITH_PETSC OR NOT DEAL_II_WITH_SLEPC OR DEAL_II_PETSC_WITH_COMPLEX) # keep in one line
  MESSAGE(FATAL_ERROR "
Error! This application requires a deal.II library that was configured with the following options:
    DEAL_II_WITH_PETSC = ON
    DEAL_II_WITH_SLEPC = ON
    DEAL_II_PETSC_WITH_COMPLEX = OFF
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options
    DEAL_II_WITH_PETSC = ${DEAL_II_WITH_PETSC}
    DEAL_II_WITH_SLEPC = ${DEAL_II_WITH_SLEPC}
    DEAL_II_PETSC_WITH_COMPLEX = ${DEAL_II_PETSC_WITH_COMPLEX}
which conflict with the requirements."
    )
ENDIF()


# Search for Boost.
# ------------------------------------------
find_package(Boost REQUIRED COMPONENTS program_options)
set(Boost_USE_STATIC_LIBS ON)


# Setup targets
# =========================================================

# Library
# ---------------------------------------------------------
# file(GLOB SOURCE_FILES "${PROJECT_SOURCE_DIR}/src/*.cc")
set(LIBRARY_SOURCE_FILES
      "models/base_model.cc"
      "models/quantum_well.cc"
      "utilities.cc"
      "parameters_parsing.cc"
      "output_results.cc"
      # "stationary_schrodinger.cc"
      "kohn_sham.cc"
      "hartree.cc"
      "dft.cc"
      "xc.cc"
      # "time_dependent.cc"
      "external_potentials.cc"
   )
# Add prefix to all items in list.
list(TRANSFORM LIBRARY_SOURCE_FILES PREPEND "lib/src/")

add_library(${LIBRARY} STATIC ${LIBRARY_SOURCE_FILES})
DEAL_II_SETUP_TARGET(${LIBRARY})
target_link_libraries(${LIBRARY} ${Boost_LIBRARIES})
target_include_directories( ${LIBRARY}
    PUBLIC
      "${PROJECT_SOURCE_DIR}/lib/include"
      "${PROJECT_SOURCE_DIR}/bundled/magic_enum"
  )

# set_target_properties(
#     ${TARGET} PROPERTIES
#     CXX_STANDARD 17
#     CXX_STANDARD_REQUIRED YES
#     CXX_EXTENSIONS NO
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
#     LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
#     ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
# )


# Executables
# ---------------------------------------------------------
foreach(target ${TARGETS})
    add_executable(${target} src/${target}.cc)
    DEAL_II_SETUP_TARGET(${target})
    target_include_directories( ${target}
        PUBLIC
          "${PROJECT_SOURCE_DIR}/lib/include"
          "${PROJECT_SOURCE_DIR}/bundled/magic_enum"
      )
    target_link_libraries( ${target}
        femdft
        ${Boost_LIBRARIES}
      )
endforeach()

# vim: ts=4 sts=4 sw=4
