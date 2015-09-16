# Overrides contents of all variables previously set by CMake
if(NOT DEFINED ENV{CXXFLAGS})
    if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
        # Compilation of Boost uncovers some bugs with Intel's support for C++11
        # For Intel compilers older that 14.0.0 continue using -std=gnu++98
        execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE ICPC_VERSION)
        if(ICPC_VERSION VERSION_LESS 14.0.0)
            set(CMAKE_CXX_FLAGS "-fPIC -std=gnu++98")
        else()
            if(HAS_CXX11_SUPPORT)
                set(CMAKE_CXX_FLAGS "-fPIC ${CXX11_COMPILER_FLAGS}")
            else()
                set(CMAKE_CXX_FLAGS "-fPIC -std=gnu++98")
            endif()
        endif()
        set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -debug -DDEBUG -Wall -Wuninitialized -Wno-unknown-pragmas")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
        # FIXME not sure this is actually needed...
        set(CMAKE_CXX_LINK_FLAGS    "-shared-intel")
    endif()
endif()
