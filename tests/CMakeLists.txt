cmake_policy(SET CMP0063 NEW)

# Prepare "Catch" library for other executables
set(CATCH_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external/Catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(unit_tests unit_tests.cpp)
set_target_properties(unit_tests
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )
target_compile_definitions(unit_tests
  PRIVATE
    "$<$<BOOL:${HAS_CXX11}>:HAS_CXX11>"
    "$<$<NOT:$<BOOL:${HAS_CXX11}>>:CATCH_CONFIG_NO_CPP11>"
  )

add_subdirectory(benchmark)
add_subdirectory(bi_operators)
add_subdirectory(cpcm)
add_subdirectory(dielectric_profile)
add_subdirectory(gepol)
add_subdirectory(green)
add_subdirectory(iefpcm)
add_subdirectory(input)
add_subdirectory(mmfq)
add_subdirectory(numerical_quadrature)
add_subdirectory(utils)

target_link_libraries(unit_tests
  PRIVATE
    pcm-static
    Catch
  )

add_subdirectory(C_host)
if(TEST_Fortran_API)
  add_subdirectory(Fortran_host)
endif()
if(BUILD_STANDALONE)
  add_subdirectory(standalone)
endif()
