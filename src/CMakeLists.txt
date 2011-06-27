add_executable ("pcmsolver.x" test.cpp)

include_directories (solverlib)
include_directories (pedralib)
include_directories (cavitylib)
include_directories (metal)
include_directories (green)

message ("CBLAS libraries in /CMakeLists.txt ${CBLAS_LIBRARIES}")
message ("LIBS in CMakeLists.txt ${LIBS}")

target_link_libraries (pcmsolver.x green metalsphere pedra cavity solver ${LIBS} ${CBLAS_LIBRARIES}) 

add_subdirectory(solverlib)
add_subdirectory(cavitylib)
add_subdirectory(pedralib)
add_subdirectory(metal)
add_subdirectory(green)
