include_directories(${PROJECT_SOURCE_DIR}/api)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

file(COPY ${PROJECT_SOURCE_DIR}/api/pcmsolver.h DESTINATION ${PROJECT_BINARY_DIR}/include)

if(ENABLE_FORTRAN_API)
endif()
