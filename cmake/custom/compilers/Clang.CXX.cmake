# Overrides contents of all variables previously set by CMake
if(NOT DEFINED ENV{CXXFLAGS})
    if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
        # Discover C++11 support
        set(CXX_STANDARD_FLAG "-std=gnu++98")
        if(ENABLE_CXX11_SUPPORT)
            include(CheckCXX11)
            discover_cxx11_support(CXX_STANDARD_FLAG)
        endif()

        set(CMAKE_CXX_FLAGS "${CXX_STANDARD_FLAG} -Qunused-arguments -fcolor-diagnostics")
        set(CMAKE_CXX_FLAGS_DEBUG    "-O0 -DDEBUG -Wall -Wextra -Winit-self -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Wwrite-strings -Weffc++ -Wdocumentation")
        set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG -Wno-unused")
    endif()
endif()
