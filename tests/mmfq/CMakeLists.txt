add_library(mmfq-tests OBJECT mmfq-ohno.cpp)
if(BUILD_CUSTOM_BOOST)
  add_dependencies(mmfq-tests custom_boost)
endif()

# Copy reference files to ${PROJECT_BINARY_DIR}/tests/mmfq (aka ${CMAKE_CURRENT_BINARY_DIR})
list(APPEND reference_files FQ-2xH2O.npy)
file(COPY ${reference_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# mmfq.cpp test
add_Catch_test(mmfq-ohno "mmfq;ohno")

