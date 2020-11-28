set(CMAKE_C_FLAGS         "${CMAKE_C_FLAGS} -std=c99 -DRESTRICT=restrict -DFUNDERSCORE=1 -Qunused-arguments -fcolor-diagnostics")
set(CMAKE_C_FLAGS_DEBUG   "-O0 -DDEBUG -g3 -Wall -Wextra -Winit-self -Wuninitialized -Wmissing-declarations -Wwrite-strings -Wno-sign-compare")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
