set(tsless_SOURCES
	tsless_cavity_interface.F
  )

add_library(tsless ${tsless_SOURCES})

install(TARGETS tsless ARCHIVE DESTINATION lib)
