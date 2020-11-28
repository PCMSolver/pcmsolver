set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -debug -DDEBUG -Wall -Wuninitialized -Wno-unknown-pragmas -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
