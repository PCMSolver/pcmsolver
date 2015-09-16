message(STATUS "System                : ${CMAKE_SYSTEM_NAME}")
message(STATUS "Processor type        : ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message(STATUS "C++ compiler flags    : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${cmake_build_type_toupper}} ${CXX_ARCHITECTURE_FLAGS}")
message(STATUS "C compiler flags      : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${cmake_build_type_toupper}} ${C_ARCHITECTURE_FLAGS}")
message(STATUS "Fortran compiler flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${cmake_build_type_toupper}} ${Fortran_ARCHITECTURE_FLAGS}")

get_directory_property(LIST_OF_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
message(STATUS "Definitions           : ${LIST_OF_DEFINITIONS}")
unset(LIST_OF_DEFINITIONS)
