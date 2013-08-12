set(green_SOURCES
	GreensFunction.cpp
	GreensFunctionFactory.cpp
	Vacuum.cpp
	UniformDielectric.cpp
	IonicLiquid.cpp
   )

add_library(green ${green_SOURCES})
add_dependencies(green libtaylor)        		

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${green_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)
                                         
install(TARGETS green ARCHIVE DESTINATION lib)		

install(FILES GreensFunction.hpp DESTINATION include)         			
install(FILES GreensFunctionFactory.hpp DESTINATION include)         			
install(FILES Vacuum.hpp DESTINATION include)                 			
install(FILES UniformDielectric.hpp DESTINATION include)      			
install(FILES IonicLiquid.hpp DESTINATION include)      			
