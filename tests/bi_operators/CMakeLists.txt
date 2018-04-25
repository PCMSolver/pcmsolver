target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/bi_operators_collocation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bi_operators_purisima.cpp
  )

# bi_operators_purisima.cpp test
add_Catch_test(
  NAME
    bi_operators_purisima
  LABELS
    bi_operators
    bi_operators_purisima
  REFERENCE_FILES
    uniformdielectric_D_purisima.npy
    vacuum_D_purisima.npy
  )

# bi_operators_collocation.cpp test
add_Catch_test(
  NAME
    bi_operators_collocation
  LABELS
    bi_operators
    bi_operators_collocation
  REFERENCE_FILES
    sphericalsharp_D_collocation.npy
    sphericalsharp_S_collocation.npy
    erfsphericaldiffuse_D_collocation.npy
    erfsphericaldiffuse_S_collocation.npy
    logsphericaldiffuse_D_collocation.npy
    logsphericaldiffuse_S_collocation.npy
    tanhsphericaldiffuse_D_collocation.npy
    tanhsphericaldiffuse_S_collocation.npy
    uniformdielectric_S_collocation.npy
    vacuum_S_collocation.npy
    uniformdielectric_D_collocation.npy
    vacuum_D_collocation.npy
  )

# This executable updates the .npy files containing the reference values
if(BUILD_STANDALONE)
  add_executable(update_reference_files update_reference_files.cpp)
  target_link_libraries(update_reference_files pcm-static)
  target_include_directories(update_reference_files
    PRIVATE
      ${PROJECT_SOURCE_DIR}/tests
    )
endif()
