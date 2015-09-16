# List of headers
list(APPEND headers_list BoundingBox.hpp Compression.hpp ConAnsatzFunction.hpp Constants.hpp ElementTree.hpp Energy.hpp GenericAnsatzFunction.hpp LinAnsatzFunction.hpp Mask.hpp Randwerte.hpp SparseMatrix.hpp Topology.hpp Transformations.hpp Volume.hpp WEM.hpp Wavelet.hpp readPoints.hpp)

# List of sources
list(APPEND sources_list ConAnsatzFunction.cpp Energy.cpp GenericAnsatzFunction.cpp LinAnsatzFunction.cpp Mask.cpp SparseMatrix.cpp Topology.cpp Transformations.cpp Volume.cpp Wavelet.cpp readPoints.cpp)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(wemFiles STATIC ${sources_list})
# Specify dependencies for the library (if any)
#add_dependencies(wemFiles )
set_property(GLOBAL APPEND PROPERTY PCMSolver_LIBRARIES wemFiles)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(wemFiles custom_boost)
endif()
install(TARGETS wemFiles ARCHIVE DESTINATION lib)

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${PROJECT_SOURCE_DIR}/wemFiles)
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/wemFiles)
endforeach()
