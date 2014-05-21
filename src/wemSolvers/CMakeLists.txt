set(wemSolvers_SOURCES
  GaussLegendre.cpp
  GaussSquare.cpp
  WEMPCG.cpp
  WEMPGMRES.cpp
  WEMRHS.cpp
)

if(ENABLE_TESTS)
  add_library(utils ${wemFiles_SOURCES})
  #add_dependencies(wemFiles libgetkw)
  install(TARGETS wemFiles ARCHIVE DESTINATION lib)
endif()

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${wemSolvers_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)

install(FILES Cubature.hpp DESTINATION include)
install(FILES GaussLegendre.hpp DESTINATION include)
install(FILES GaussSquare.hpp DESTINATION include)
install(FILES Quadrature.hpp DESTINATION include)
install(FILES WEMPCG.hpp DESTINATION include)
install(FILES WEMPGMRES.hpp DESTINATION include)
install(FILES WEMRHS.hpp DESTINATION include)
