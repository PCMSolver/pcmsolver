#.rst:
#
# Enables extended diagnostics flags for the GNU C++ compiler
#
# Variables used::
#
#   ENABLE_EXTENDED_DIAGNOSTICS
#
# Variables modified::
#
#   EXDIAG_CXX_FLAGS
#
# autocmake.cfg configuration::
#
#   docopt: --exdiag Enable C++ extended diagnostics flags [default: False].
#   define: '-DENABLE_EXTENDED_DIAGNOSTICS=%s' % arguments['--exdiag']

option(ENABLE_EXTENDED_DIAGNOSTICS "Enable extendend diagnostics compiler flags" OFF)

if(ENABLE_EXTENDED_DIAGNOSTICS)
  set(EXDIAG_CXX_FLAGS "-Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wsuggest-final-types -Wsuggest-final-methods -Wsuggest-override -Wuseless-cast -Wunsafe-loop-optimizations")
endif()
