option(BUILD_C_TEST "Build the C test program" OFF)

add_library(cgetkw 
	getkw.c
	pygetkw.c
	)

install(TARGETS cgetkw ARCHIVE DESTINATION lib)
install(FILES getkw.h 
    DESTINATION include
	)

if (BUILD_C_TEST)
	add_executable(getkw_test getkw_test.c)
	target_link_libraries(getkw_test cgetkw)
endif (BUILD_C_TEST)
