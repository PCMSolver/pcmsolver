include_directories(${CMAKE_CURRENT_LIST_DIR}/cavity
                    ${CMAKE_CURRENT_LIST_DIR}/green
                    ${CMAKE_CURRENT_LIST_DIR}/green/dielectric_profile
                    ${CMAKE_CURRENT_LIST_DIR}/interface
                    ${CMAKE_CURRENT_LIST_DIR}/bi_operators
                    ${CMAKE_CURRENT_LIST_DIR}/solver
                    ${CMAKE_CURRENT_LIST_DIR}/utils
                   )

add_subdirectory(cavity)
add_subdirectory(green)
add_subdirectory(interface)
add_subdirectory(metal)
add_subdirectory(bi_operators)
add_subdirectory(pedra)
add_subdirectory(solver)
add_subdirectory(utils)

get_property(_pcmsolver_cxx_sources GLOBAL PROPERTY PCMSolver_CXX_SOURCES)
get_property(_pcmsolver_c_sources GLOBAL PROPERTY PCMSolver_C_SOURCES)
get_property(_pcmsolver_fortran_sources GLOBAL PROPERTY PCMSolver_Fortran_SOURCES)
add_library(pcm ${_pcmsolver_cxx_sources} ${_pcmsolver_c_sources} ${_pcmsolver_fortran_sources})
add_dependencies(pcm libgetkw)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(pcm libgetkw)
   add_dependencies(pcm custom_boost)
endif()
install(TARGETS pcm ARCHIVE DESTINATION lib)

# This directory contains source file for executables
if(BUILD_STANDALONE)
   add_subdirectory(bin)
endif()
