# List of headers
list(APPEND headers_list BoundingBox.hpp Compression.hpp ConAnsatzFunction.hpp Constants.hpp ElementTree.hpp Energy.hpp GenericAnsatzFunction.hpp LinAnsatzFunction.hpp Mask.hpp Randwerte.hpp SparseMatrix.hpp Topology.hpp Transformations.hpp Volume.hpp WEM.hpp Wavelet.hpp readPoints.hpp)

# List of sources
list(APPEND sources_list ConAnsatzFunction.cpp Energy.cpp GenericAnsatzFunction.cpp LinAnsatzFunction.cpp Mask.cpp SparseMatrix.cpp Topology.cpp Transformations.cpp Volume.cpp Wavelet.cpp readPoints.cpp)

# Update bar chart
if(BUILD_DOCUMENTATION)
    update_bar_chart(${CMAKE_CURRENT_LIST_DIR})
endif()

set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(_source ${sources_list})
    set_property(GLOBAL APPEND PROPERTY PCMSolver_CXX_SOURCES ${CMAKE_CURRENT_LIST_DIR}/${_source})
endforeach()
# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include/wemFiles)
endforeach()
