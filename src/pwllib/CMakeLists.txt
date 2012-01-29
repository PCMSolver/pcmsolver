set (pwl_C_SOURCES 
  WEMPCG_pwl.c
  WEMPGMRES_pwl.c
  WEMRHS_pwl.c
  WEM_pwl.c
  basis_pwl.c
  compression_pwl.c
  det.c
  dwt_pwl.c
  energy_pwl.c
  integrate_pwl.c
  interpolate_pwl.c
  intlin1.c
  intlin2.c
  intlin3.c
  intlin4.c
  intvector_pwl.c
  main_pwl.c
  mask_pwl.c
  phi.c
  postproc_pwl.c
  precond_pwl.c
  topology_pwl.c
  trafos_pwl.c
  volume_pwl.c
  )

add_library(pwl ${pwl_C_SOURCES} ${pwl_CXX_SOURCES})

install(FILES WEMPCG_pwl.h       DESTINATION include)
install(FILES WEMPGMRES_pwl.h    DESTINATION include)
install(FILES WEMRHS_pwl.h       DESTINATION include)
install(FILES WEM_pwl.h          DESTINATION include)
install(FILES basis_pwl.h        DESTINATION include)
install(FILES compression_pwl.h  DESTINATION include)
install(FILES dwt_pwl.h          DESTINATION include)
install(FILES energy_pwl.h       DESTINATION include)
install(FILES integrate_pwl.h    DESTINATION include)
install(FILES interpolate_pwl.h  DESTINATION include)
install(FILES intlin1.h          DESTINATION include)
install(FILES intlin2.h          DESTINATION include)
install(FILES intlin3.h          DESTINATION include)
install(FILES intlin4.h          DESTINATION include)
install(FILES intvector_pwl.h    DESTINATION include)
install(FILES mask_pwl.h         DESTINATION include)
install(FILES phi.h              DESTINATION include)
install(FILES postproc_pwl.h     DESTINATION include)
install(FILES precond_pwl.h      DESTINATION include)
install(FILES topology_pwl.h     DESTINATION include)
install(FILES trafos_pwl.h       DESTINATION include)
install(FILES volume_pwl.h       DESTINATION include)

install(TARGETS pwl
    ARCHIVE DESTINATION lib)

