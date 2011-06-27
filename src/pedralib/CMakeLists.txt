set(pedra_SOURCES 
    mywlkdin.f90
    gp_dblas1.F            
    makecav.F              
    printpkg.F             
    sircav.F               
    sirpcmder.F            
    subroutines_cav.F      
    traditionalPCMcavity.F
    pcmcav.F
)

add_library(pedra ${pedra_SOURCES})

install(TARGETS pedra
	ARCHIVE DESTINATION lib)
