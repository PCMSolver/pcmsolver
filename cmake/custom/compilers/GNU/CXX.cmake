set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3 -DDEBUG -Wall -Wextra -Winit-self -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Wwrite-strings -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wno-unused")
