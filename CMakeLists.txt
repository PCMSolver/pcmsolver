cmake_minimum_required (VERSION 2.6)
project (pcmsolver)

# build executable
set (EXECUTABLE "pcmsolver.exe")


set (pcmsolver_SOURCES
pcmsolver.cpp
TraditionalPCMOperator.cpp 
TraditionalPCMOperator.h)


add_definitions(-I/home/notasoul/luca/eigen)

add_executable ("pcmsolver.exe" ${pcmsolver_SOURCES})

# install executables and scripts
install (TARGETS ${EXECUTABLE} 
         RUNTIME DESTINATION "bin")

install (PROGRAMS ${SCRIPTS}
         DESTINATION "bin") 



