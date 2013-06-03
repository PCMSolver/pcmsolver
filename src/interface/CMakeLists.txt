set(interface_SOURCES 
	#	CavityFactory.cpp
	#	GreensFunctionFactory.cpp
	#	SolverFactory.cpp
        Interface.cpp
)

add_library(interface ${interface_SOURCES})
add_dependencies(interface external)

install(TARGETS interface 
	ARCHIVE DESTINATION lib)

#install(FILES CavityFactory.h DESTINATION include)
#install(FILES GreensFunctionFactory.h DESTINATION include)
#install(FILES SolverFactory.h DESTINATION include)
install(FILES Interface.h DESTINATION include)
