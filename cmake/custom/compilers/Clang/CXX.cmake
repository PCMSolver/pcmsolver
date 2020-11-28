set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments -fcolor-diagnostics")
set(CMAKE_CXX_FLAGS_DEBUG    "-O0 -DDEBUG -Wall -Wextra -Winit-self -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Wwrite-strings -Weffc++ -Wdocumentation -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG -Wno-unused")
