set(green_SOURCES
	GreensFunction.cpp
	GreensFunctionFactory.cpp
	Vacuum.cpp
	UniformDielectric.cpp
	IonicLiquid.cpp
   )

add_library(green ${green_SOURCES})
add_dependencies(green libtaylor)        		
                                         
install(TARGETS green ARCHIVE DESTINATION lib)		

install(FILES GreensFunction.hpp DESTINATION include)         			
install(FILES GreensFunctionFactory.hpp DESTINATION include)         			
install(FILES Vacuum.hpp DESTINATION include)                 			
install(FILES UniformDielectric.hpp DESTINATION include)      			
install(FILES IonicLiquid.hpp DESTINATION include)      			
