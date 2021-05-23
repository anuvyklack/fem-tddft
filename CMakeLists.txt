# PROJECT(schrodinger VERSION 0.0.3)
PROJECT(fem-tddft)

SET(TARGET tddft)

# SET(DEAL_II_DIR ~/software/dealii-9.2.0/build)


# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Search for deal.II library.
# ---------------------------
FIND_PACKAGE(deal.II 9.2.0 QUIET
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


# Search for Boost::program_options library.
# ------------------------------------------
find_package(Boost REQUIRED COMPONENTS program_options)
set(Boost_USE_STATIC_LIBS ON)


# Setup targets
# -------------
# file(GLOB SOURCE_FILES "${PROJECT_SOURCE_DIR}/src/*.cc")
set(SOURCE_FILES
      "src/main.cc"
      "src/utilities.cc"
      "src/parameters_parsing.cc"
      "src/model.cc"
      "src/output_results.cc"
      # "src/stationary_schrodinger.cc"
      "src/kohn_sham.cc"
      "src/hartree.cc"
      "src/dft.cc"
      # "src/time_dependent.cc"
   )
add_executable(${TARGET} ${SOURCE_FILES})
target_link_libraries(${TARGET} ${Boost_LIBRARIES})
target_include_directories(
    ${TARGET} PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/bundled/magic_enum"
  )

# set_target_properties(${TARGET} PROPERTIES
#     CXX_STANDARD 17
#     CXX_STANDARD_REQUIRED YES
#     CXX_EXTENSIONS NO
# )

DEAL_II_SETUP_TARGET(${TARGET})


# # Test target
# # -----------------
# set(TEST_SOURCE_FILES "src/test.cc"
#                       "src/general.cc"
#                       "src/model.cc"
#                       "src/ground_state.cc"
#    )
# add_executable(test ${TEST_SOURCE_FILES})
# target_link_libraries(test ${Boost_LIBRARIES})
# target_include_directories(test PUBLIC
#                            "${PROJECT_SOURCE_DIR}/include")
# DEAL_II_SETUP_TARGET(test)


# vim: ts=2 sts=2 sw=2
