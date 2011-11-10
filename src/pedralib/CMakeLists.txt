set(pedra_SOURCES 
  pcm_utils.F
  pcm_print.F
  pcm_inertia.F90
  pcm_der.F
  pcm_cav.F
  pcm_cav_interface.F
  makecav.F
  pcm_dblas1.F
  )

add_library(pedra ${pedra_SOURCES})

install(TARGETS pedra
  ARCHIVE DESTINATION lib)
