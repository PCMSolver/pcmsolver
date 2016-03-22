add_library(bi_operators-tests OBJECT bi_operators_collocation.cpp bi_operators_numerical.cpp bi_operators_purisima.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(bi_operators-tests custom_boost)
endif()

# Copy reference files to ${PROJECT_BINARY_DIR}/tests/bi_operators (aka ${CMAKE_CURRENT_BINARY_DIR})
list(APPEND reference_files tanhsphericaldiffuse_D_collocation.npy uniformdielectric_D_purisima.npy vacuum_D_purisima.npy
                            tanhsphericaldiffuse_S_collocation.npy uniformdielectric_S_collocation.npy vacuum_S_collocation.npy
                            uniformdielectric_D_collocation.npy vacuum_D_collocation.npy)
file(COPY ${reference_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# bi_operators_purisima.cpp test
add_Catch_test(bi_operators_purisima "bi_operators;bi_operators_purisima")

# bi_operators_numerical.cpp test
add_Catch_test(bi_operators_numerical "bi_operators;bi_operators_numerical")

# bi_operators_collocation.cpp test
add_Catch_test(bi_operators_collocation "bi_operators;bi_operators_collocation")

# This executable updates the .npy files containing the reference values
if(BUILD_STANDALONE)
  link_directories(${PROJECT_BINARY_DIR}/lib)
  add_executable(update_reference_files update_reference_files.cpp)
  if(BUILD_CUSTOM_BOOST)
    add_dependencies(update_reference_files custom_boost)
  endif()
  if(STATIC_LIBRARY_ONLY)
    target_link_libraries(update_reference_files pcm-static ${ZLIB_LIBRARIES})
  else()
    target_link_libraries(update_reference_files pcm-shared)
  endif()
endif()
