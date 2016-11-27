option(ENABLE_CXX11_SUPPORT "Enable C++11 compiler support" ON)

set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

include(GNU.CXX)
include(Clang.CXX)
include(Intel.CXX)
