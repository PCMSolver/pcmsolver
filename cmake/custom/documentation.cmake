include(FindDoxygen QUIET)

if(DOXYGEN_FOUND)
    configure_file(doc/Doxyfile.in Doxyfile @ONLY)
    # Really ugly...
    add_custom_target(
        bar_charts
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_BINARY_DIR}/bin/counter.py"
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/cavity.py"       WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/green.py"        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/interface.py"    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/metal.py"        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/bi_operators.py" WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/pedra.py"        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/solver.py"       WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/total.py"        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/utils.py"        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/gfx
    )
    add_custom_target(doc
        COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_BINARY_DIR}/bin/counter.py"
        COMMAND ${DOXYGEN_EXECUTABLE} WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
    add_dependencies(doc bar_charts)
#   add_custom_target(update_gh-pages
#       COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_BINARY_DIR}/bin/update_gh-pages.py"
#       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
#   )
#   add_dependencies(update_gh-pages doc)
else()
	message(STATUS "Doxygen missing. You won't be able to create docs.")
endif()
