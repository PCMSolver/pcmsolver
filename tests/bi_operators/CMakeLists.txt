add_library(bi_operators-tests OBJECT
  bi_operators_collocation.cpp
  bi_operators_purisima.cpp
  )
if(BUILD_CUSTOM_BOOST)
  add_dependencies(bi_operators-tests custom_boost)
endif()

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
  if(BUILD_CUSTOM_BOOST)
    add_dependencies(update_reference_files custom_boost)
  endif()
endif()
