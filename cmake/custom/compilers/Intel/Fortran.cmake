set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -w -fpp -assume byterecl -traceback -nosave")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -warn all")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
