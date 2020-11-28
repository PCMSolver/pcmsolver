set(CMAKE_C_FLAGS         "${CMAKE_C_FLAGS} -std=c99 -DRESTRICT=restrict -DFUNDERSCORE=1 -fPIC")
set(CMAKE_C_FLAGS_DEBUG   "-O0 -g3 -DDEBUG -Wall -Wextra -Winit-self -Wuninitialized -Wmissing-declarations -Wwrite-strings -Wno-sign-compare")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
