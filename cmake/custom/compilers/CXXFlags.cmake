option_with_print(ENABLE_CXX11_SUPPORT "Enable C++11 compiler support" ON)

# Discover C++11 support
if(ENABLE_CXX11_SUPPORT)
  include(SetCompilerFlag)
  set_compiler_flag(
    RESULT cxx11_flag
    LANGUAGE CXX
    #REQUIRED
    FLAGS "-std=c++11;/std:c++11;-std=c++0x;--c++11;--c++0x"
    )
  if(cxx11_flag)
    set(CXX_STANDARD_FLAG ${cxx11_flag})
    set(HAS_CXX11 TRUE)
    message(STATUS "Using C++11 standard")
  else()
    set(CXX_STANDARD_FLAG "-std=gnu++98")
    set(HAS_CXX11 FALSE)
    message(STATUS "Using C++03 standard")
    message(WARNING "Upgrade your compiler! C++03 support is DEPRECATED and will be removed")
  endif()
else()
  message(WARNING "Upgrade your compiler! C++03 support is DEPRECATED and will be removed")
endif()

include(GNU.CXX)
include(Clang.CXX)
include(Intel.CXX)
