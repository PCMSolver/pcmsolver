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

install(FILES GreensFunction.h DESTINATION include)         			
install(FILES GreensFunctionFactory.h DESTINATION include)         			
install(FILES Vacuum.h DESTINATION include)                 			
install(FILES UniformDielectric.h DESTINATION include)      			
install(FILES IonicLiquid.h DESTINATION include)      			
