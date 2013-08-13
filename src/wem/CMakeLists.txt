set (wem_SOURCES 
  WEM.c
  WEMPCG.c
  WEMPGMRES.c
  WEMRHS.c
  basis.c
  compression.c
  constants.c
  dwt.c
  energy.c
  gauss_legendre.c
  gauss_square.c
  integrate.c
  interpolate2.c
  intkon1.c
  intkon2.c
  intkon3.c
  intkon4.c
  intvector.c
  kern.c
  mask.c
  molecule.c
  postproc.c
  precond.c
  randwerte.c
  read_points_asc.c
  sparse.c
  sparse2.c
  topology.c
  trafos.c
  vector2.c
  vector3.c
  volume.c
  )

add_library(wem ${wem_SOURCES})

# Export the sources up to src directory level
set(local_sources "")
foreach(source ${wem_SOURCES})
	get_filename_component(sourcepath ${source} ABSOLUTE)
  	set(local_sources ${sourcepath} ${local_sources})
endforeach()
set(SOURCES ${local_sources} ${SOURCES} PARENT_SCOPE)

install(TARGETS wem ARCHIVE DESTINATION lib)

install(FILES WEM.h            DESTINATION include)
install(FILES WEMPCG.h         DESTINATION include)
install(FILES WEMPGMRES.h      DESTINATION include)
install(FILES WEMRHS.h         DESTINATION include)
install(FILES basis.h          DESTINATION include)
install(FILES compression.h    DESTINATION include)
install(FILES constants.h      DESTINATION include)
install(FILES cubature.h       DESTINATION include)
install(FILES dwt.h            DESTINATION include)
install(FILES energy.h         DESTINATION include)
install(FILES gamma.h          DESTINATION include)
install(FILES gauss_legendre.h DESTINATION include)
install(FILES gauss_square.h   DESTINATION include)
install(FILES init_points.h    DESTINATION include)
install(FILES integrate.h      DESTINATION include)
install(FILES interpolate.h    DESTINATION include)
install(FILES intkon1.h        DESTINATION include)
install(FILES intkon2.h        DESTINATION include)
install(FILES intkon3.h        DESTINATION include)
install(FILES intkon4.h        DESTINATION include)
install(FILES intvector.h      DESTINATION include)
install(FILES kern.h           DESTINATION include)
install(FILES mask.h           DESTINATION include)
install(FILES molecule.h       DESTINATION include)
install(FILES postproc.h       DESTINATION include)
install(FILES pot.h            DESTINATION include)
install(FILES precond.h        DESTINATION include)
install(FILES quadrature.h     DESTINATION include)
install(FILES randwerte.h      DESTINATION include)
install(FILES read_points.h    DESTINATION include)
install(FILES sparse.h         DESTINATION include)
install(FILES sparse2.h        DESTINATION include)
install(FILES topology.h       DESTINATION include)
install(FILES trafos.h         DESTINATION include)
install(FILES vector2.h        DESTINATION include)
install(FILES vector3.h        DESTINATION include)
install(FILES volume.h         DESTINATION include)
