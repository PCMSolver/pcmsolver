set(Doxygen_FIND_QUIETLY TRUE)
include(FindDoxygen)

if(DOXYGEN_FOUND)
    configure_file(
        ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
        ${PROJECT_BINARY_DIR}/Doxyfile)
    # Really ugly...
    add_custom_target(
	bar_charts
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/cavity.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/green.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/interface.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/metal.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/operators_diagonal.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/pedra.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/pwl.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/solver.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/total.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/tsless.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/utils.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/wavcav.py" 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/doc/gfx/matplotlib/wem.py") 
    add_custom_target(
        doc 
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_BINARY_DIR}/bin/counter.py"
        COMMAND ${DOXYGEN_EXECUTABLE} WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
    add_dependencies(doc bar_charts)
    add_custom_target(
	update_gh-pages
	COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_BINARY_DIR}/bin/update_gh-pages.py"
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
    add_dependencies(update_gh-pages doc)
else()
	message(STATUS "Doxygen missing. You won't be able to create docs.")
endif()
