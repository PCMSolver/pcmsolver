set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fimplicit-none -fautomatic -fmax-errors=5")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbacktrace -Wall -fcheck=all")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops -ftree-vectorize")
