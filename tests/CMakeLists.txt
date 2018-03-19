cmake_policy(SET CMP0063 NEW)

include_directories(
  SYSTEM
    ${PROJECT_SOURCE_DIR}/external/Catch
  )
include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/tests
  )

add_subdirectory(input)
add_subdirectory(numerical_quadrature)
add_subdirectory(gepol)
add_subdirectory(dielectric_profile)
add_subdirectory(green)
add_subdirectory(bi_operators)
add_subdirectory(cpcm)
add_subdirectory(iefpcm)
add_subdirectory(utils)
add_subdirectory(benchmark)

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
target_compile_definitions(unit_tests
  PRIVATE
    $<$<NOT:$<BOOL:${ENABLE_CXX11_SUPPORT}>>:CATCH_CONFIG_NO_CPP11>
  )
target_link_libraries(unit_tests
  PRIVATE
    pcm-static
  )

add_subdirectory(C_host)
if(ENABLE_Fortran_API)
  add_subdirectory(Fortran_host)
endif()
if(BUILD_STANDALONE)
  add_subdirectory(standalone)
endif()
