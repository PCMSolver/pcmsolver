set(green_source
    GreensFunctionInterface.cpp
    GreensFunction.cpp
    Vacuum.cpp
    MetalSphere.cpp
    UniformDielectric.cpp
    GreensFunctionSum.cpp
)

add_library(green ${green_source})

install(TARGETS green 
	ARCHIVE DESTINATION lib)

install(FILES GreensFunctionInterface.h DESTINATION include)
install(FILES GreensFunction.h DESTINATION include)
install(FILES MetalSphere.h DESTINATION include)
install(FILES Vacuum.h DESTINATION include)
install(FILES PlanarInterface.h DESTINATION include)
install(FILES UniformDielectric.h DESTINATION include)
install(FILES GreensFunctionSum.h DESTINATION include)
