set(CMAKE_C_FLAGS         "${CMAKE_C_FLAGS} -restrict -DRESTRICT=restrict -std=c99 -fPIC")
set(CMAKE_C_FLAGS_DEBUG   "-O0 -DDEBUG -g -w3 -Wall -Wuninitialized -Wno-sign-compare")
set(CMAKE_C_FLAGS_RELEASE "-O3 -ip -DNDEBUG")
