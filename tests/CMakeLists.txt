include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_LIST_DIR})

if(NOT ENABLE_CXX11_SUPPORT)
  add_definitions(-DCATCH_CONFIG_NO_CPP11)
endif()

add_subdirectory(input)
add_subdirectory(numerical_quadrature)
add_subdirectory(gepol)
add_subdirectory(dielectric_profile)
add_subdirectory(green)
add_subdirectory(bi_operators)
add_subdirectory(cpcm)
add_subdirectory(iefpcm)
add_subdirectory(utils)

link_directories(${PROJECT_BINARY_DIR}/lib)
add_executable(unit_tests unit_tests.cpp
                          $<TARGET_OBJECTS:bi_operators-tests>
                          $<TARGET_OBJECTS:cpcm-tests>
                          $<TARGET_OBJECTS:dielectric_profile-tests>
                          $<TARGET_OBJECTS:gepol-tests>
                          $<TARGET_OBJECTS:green-tests>
                          $<TARGET_OBJECTS:iefpcm-tests>
                          $<TARGET_OBJECTS:input-tests>
                          $<TARGET_OBJECTS:numerical_quadrature-tests>
                          $<TARGET_OBJECTS:utils-tests>
                          )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(unit_tests custom_boost)
endif()
if(STATIC_LIBRARY_ONLY)
  target_link_libraries(unit_tests pcm-static)
else()
  target_link_libraries(unit_tests pcm-shared)
endif()

add_subdirectory(C_host)
if(ENABLE_FORTRAN_API)
  add_subdirectory(Fortran_host)
endif()
if(BUILD_STANDALONE)
  add_subdirectory(standalone)
endif()
