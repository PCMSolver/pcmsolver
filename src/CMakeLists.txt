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
# This directory contains source file for executables
if(BUILD_STANDALONE)
   add_subdirectory(bin)
endif()
