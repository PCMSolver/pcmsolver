find_package(Boost)

include_directories(${Boost_INCLUDE_DIRS})

add_library(getkw
	Getkw.cpp
	GetkwError.cpp
	Section.cpp
	messages.cpp
	)
set_target_properties(getkw PROPERTIES POSITION_INDEPENDENT_CODE 1)

install(TARGETS getkw
	ARCHIVE DESTINATION lib)

install(FILES Getkw.h Section.h Keyword.h GetkwError.h
	DESTINATION include
	)
