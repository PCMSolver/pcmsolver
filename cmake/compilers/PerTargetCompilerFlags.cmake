macro(set_wavcav_compiler_flags)
  set(CMAKE_C_FLAGS         " ")
  set(CMAKE_C_FLAGS_DEBUG   " ")
  set(CMAKE_C_FLAGS_RELEASE " ")

  message("empty? "       "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
  # Now set the flags as you like here
  if(CMAKE_C_COMPILER_ID MATCHES GNU)
	  set(CMAKE_C_FLAGS         "-w -std=c99 -DRESTRICT=restrict -DFUNDERSCORE=1 -fPIC")
	  set(CMAKE_C_FLAGS_DEBUG   "-O0 -g3"                                              )
	  set(CMAKE_C_FLAGS_RELEASE "-O3"                                                  )
  endif()
  
  if(CMAKE_C_COMPILER_ID MATCHES Intel)
          set(CMAKE_C_FLAGS         "-w -restrict -DRESTRICT=restrict -vec-report0 -std=c99 -fPIC")
          set(CMAKE_C_FLAGS_DEBUG   "-O0 -g -vec-report"                                          )
          set(CMAKE_C_FLAGS_RELEASE "-O3 -ip"                                                     )
  endif()
  
  if(CMAKE_C_COMPILER_ID MATCHES Clang)
	  set(CMAKE_C_FLAGS         "-w -std=c99 -DRESTRICT=restrict -DFUNDERSCORE=1 -fPIC")
	  set(CMAKE_C_FLAGS_DEBUG   "-O0 -g3"                                              )
	  set(CMAKE_C_FLAGS_RELEASE "-O3"                                                  )
  endif()
  
  if(CMAKE_C_COMPILER_ID MATCHES PGI)
	  set(CMAKE_C_FLAGS         " "                                                       )
	  set(CMAKE_C_FLAGS_DEBUG   "-g -O0"                                                  )
	  set(CMAKE_C_FLAGS_RELEASE "-O3 -fast -Munroll -Mvect=idiom -c9x -DRESTRICT=restrict")
  endif()
  
  if(CMAKE_C_COMPILER_ID MATCHES XL)
	  set(CMAKE_C_FLAGS         "-qcpluscmt")
	  set(CMAKE_C_FLAGS_DEBUG   " "         )
	  set(CMAKE_C_FLAGS_RELEASE "-O3"       )
  endif()
  
  message("in set_wavcav_compiler_flags: "       "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")

endmacro(set_wavcav_compiler_flags)
