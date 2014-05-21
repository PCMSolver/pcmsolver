set(wemFiles_SOURCES
  ConAnsatzFunction.cpp
  GenericAnsatzFunction.cpp
  LinAnsatzFunction.cpp
  Mask.cpp
  Topology.cpp
  SparseMatrix.cpp
  Transformations.cpp
  Wavelet.cpp
)

if(ENABLE_TESTS)
  add_library(utils ${wemFiles_SOURCES})
  #add_dependencies(wemFiles libgetkw)
  install(TARGETS wemFiles ARCHIVE DESTINATION lib)
endif()

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${wemFiles_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)

install(FILES BoundingBox.hpp DESTINATION include)
install(FILES ConAnsatzFunction.hpp DESTINATION include)
install(FILES Constants.hpp DESTINATION include)
install(FILES ElementTree.hpp DESTINATION include)
install(FILES GenericAnsatzFunction.hpp DESTINATION include)
install(FILES LinAnsatzFunction.hpp DESTINATION include)
install(FILES Mask.hpp DESTINATION include)
install(FILES Randwerte.hpp DESTINATION include)
install(FILES SparseMatrix.hpp DESTINATION include)
install(FILES Topology.hpp DESTINATION include)
install(FILES Transformations.hpp DESTINATION include)
install(FILES Wavelet.hpp DESTINATION include)
