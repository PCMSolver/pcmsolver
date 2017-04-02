#.rst:
#
# Detect Eigen3.
# By default, use Eigen 3.3.0 as bundled with PCMSolver.
# Look in a specific search directory, if given. If nothing is found
# there, falls back to Eigen 3.3.0 bundled with PCMSolver.
#
# autocmake.yml configuration::
#
#   docopt: "--eigen=<EIGEN3_ROOT> Root directory for Eigen3 [default: '']."
#   define: "'-DEIGEN3_ROOT=\"{0}\"'.format(arguments['--eigen'])"

if(EIGEN3_ROOT)
  set(EIGEN3_INCLUDE_DIR ${EIGEN3_ROOT}/include)
  find_package(Eigen3 3.3.0)
  message(STATUS "Eigen " ${EIGEN3_VERSION} " is located here: " ${EIGEN3_INCLUDE_DIR})
  if(NOT EIGEN3_FOUND)
    set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external/eigen3/include/eigen3)
    message(STATUS "Eigen 3.3.2 is located here: " ${EIGEN3_INCLUDE_DIR})
  endif()
else()
  set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external/eigen3/include/eigen3)
  message(STATUS "Eigen 3.3.2 is located here: " ${EIGEN3_INCLUDE_DIR})
endif()
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
