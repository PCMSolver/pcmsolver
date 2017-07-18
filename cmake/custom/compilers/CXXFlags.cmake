# We require C++11 support from the compiler and standard library.
# We do not use C++ regexes so the minimal GCC version is 4.8
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNU.CXX)
include(Clang.CXX)
include(Intel.CXX)
