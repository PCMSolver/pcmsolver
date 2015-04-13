add_library(fgetkw 
	getkw.f90
	kinds.f90
	teletype.f90
	)

install(TARGETS fgetkw ARCHIVE DESTINATION lib)
