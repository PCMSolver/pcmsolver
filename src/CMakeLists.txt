# Make sure the compiler can find include files from our library.
include_directories (solverlib)

# Make sure the linker can find the Hello library once it is built.
#link_directories (/Users/luca/programs/green/src/greenlib)
#link_directories (/Users/luca/programs/green/src/metalspherelib)

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (pcmsolver.x test.cpp)

# Link the executable to the Hello library.
target_link_libraries (pcmsolver.x green metalsphere pedra solver) 

add_subdirectory(solverlib)



