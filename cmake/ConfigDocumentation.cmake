find_package(Sphinx QUIET)

if(SPHINX_FOUND)
    add_custom_target(
        html
#        COMMAND ${CMAKE_SOURCE_DIR}/doc/preprocess_sphinx ${CMAKE_SOURCE_DIR}/doc ${PROJECT_BINARY_DIR};
                sphinx-build -b html -d _build/doctrees ${CMAKE_SOURCE_DIR}/doc _build/html
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
else()
    add_custom_target(
        html
        COMMAND echo error: please install python-sphinx and python-matplotlib first
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
endif()

set(Doxygen_FIND_QUIETLY TRUE)
include(FindDoxygen)
if(DOXYGEN_FOUND)
    configure_file(
        ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
        ${PROJECT_BINARY_DIR}/Doxyfile
    )
    add_custom_target(
        doxygen
        COMMAND ${DOXYGEN_EXECUTABLE}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
endif()

if(DOXYGEN_FOUND AND SPHINX_FOUND)
    	configure_file(
        	${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
        	${PROJECT_BINARY_DIR}/Doxyfile
    	)
	set(DOXYGEN_XML_DIR ${PROJECT_BINARY_DIR}/doc/_doxygen/xml)
	configure_file(
		${CMAKE_SOURCE_DIR}/doc/conf.py.in
		${CMAKE_SOURCE_DIR}/doc/conf.py
	)  
	add_custom_target(
		docs
		COMMAND sphinx-build -b html -d doc/doctrees ${CMAKE_SOURCE_DIR}/doc doc/html # This will generate docs with Sphinx
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	)
	add_dependencies(docs doxygen)
else()
	message("Doxygen, Sphinx or both are missing...")
endif()
