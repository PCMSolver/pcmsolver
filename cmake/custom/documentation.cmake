# This script checks the availability of the tools
# needed to build the HTML documentation pages:
# - Doxygen
# - Sphinx
# - Perl (to count the lines of code)
# - PyYAML module (to postprocess the lines of code counting)
# - Matplotlib module (to generate the bar charts)
# - Breathe module (to bridge Doxygen and Sphinx)

option(ENABLE_DOCS "Enable building of documentation" ON)

if(ENABLE_DOCS)
   find_package(Doxygen)
   find_package(Sphinx)
   find_package(Perl)

   include(find_python_module)
   find_python_module(yaml)
   find_python_module(breathe)
   find_python_module(matplotlib)
endif()

set(BUILD_DOCS ON)
if(NOT ENABLE_DOCS OR NOT DOXYGEN_FOUND OR NOT SPHINX_FOUND OR NOT PERL_FOUND OR NOT PY_YAML_FOUND OR NOT PY_BREATHE_FOUND OR NOT PY_MATPLOTLIB_FOUND)
    set(BUILD_DOCS OFF)
endif()

if(BUILD_DOCS)
    message(STATUS "Added doc target to build documentation")
    add_custom_target(doc
        COMMAND sphinx-build -b html -d doc/doctrees ${PROJECT_SOURCE_DIR}/doc doc/html
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
endif()
