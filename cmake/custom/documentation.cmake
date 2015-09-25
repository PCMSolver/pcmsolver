option(BUILD_DOCUMENTATION "Build API documentation using Doxygen" ON)

set(_build_docs OFF)
if(BUILD_DOCUMENTATION)
    find_package(Doxygen QUIET)
    find_package(Perl QUIET)
    if(DOXYGEN_FOUND AND PERL_FOUND)
        set(_build_docs ON)
    endif()
endif()

if(_build_docs)
    message(STATUS "Doxygen and Perl available to build docs")
    # Configure the cloc_tools module
    configure_file(${PROJECT_SOURCE_DIR}/tools/cloc_tools.py.in ${PROJECT_BINARY_DIR}/bin/cloc_tools.py @ONLY)

    # Configure the update_gh-pages utility script
    configure_file(${PROJECT_SOURCE_DIR}/tools/update_gh-pages.py.in ${PROJECT_BINARY_DIR}/bin/update_gh-pages.py @ONLY)

    # Configure Doxyfile
    configure_file(doc/Doxyfile.in Doxyfile @ONLY)

    # Build documenation from Doxygen strings in source files
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    # Updates website
    add_custom_target(update_gh-pages
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/bin/update_gh-pages.py
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    add_dependencies(update_gh-pages doc)
endif(_build_docs)

# Updates bar charts for given directory of source files
# _src_dir is the source directory
# _lang is the language: CXX, C or Fortran
macro(update_bar_chart _src_dir _lang)
    get_filename_component(_script_name ${_src_dir} NAME)
    set(_working_dir ${PROJECT_BINARY_DIR}/doc/gfx/matplotlib)
    file(MAKE_DIRECTORY ${_working_dir})
    # Generate bar chart script
    set(_counter "import sys; \
                  sys.path.append('${PROJECT_BINARY_DIR}/bin'); \
                  from cloc_tools import bar_chart; \
                  bar_chart('${_src_dir}', '${_lang}')"
                  )
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "${_counter}"
                    WORKING_DIRECTORY ${_working_dir}
                    RESULT_VARIABLE _script_generated
                    )
    # Update bar chart
    execute_process(COMMAND ${PYTHON_EXECUTABLE} ${_working_dir}/${_script_name}.py
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx
                    RESULT_VARIABLE _chart_generated
                    )
endmacro()
