#.rst:
#
# Enable profiling with gperftools.
#
# Variables used::
#
#   ENABLE_GPERFTOOLS
#
# autocmake.yml configuration::
#
#   docopt: "--gperf Enable profiling with gperftools [default: False]."
#   define: "'-DENABLE_GPERFTOOLS={0}'.format(arguments['--gperf'])"

option_with_print(ENABLE_GPERFTOOLS "Enable profiling with gperftools" OFF)

if(ENABLE_GPERFTOOLS)
  message(STATUS "Linking against gperftools libraries for profiling")
  find_package(Gperftools COMPONENTS tcmalloc profiler)
endif()
