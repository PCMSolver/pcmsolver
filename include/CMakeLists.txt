list(APPEND headers_list Citation.hpp Cxx11Workarounds.hpp ErrorHandling.hpp LoggerInterface.hpp TimerInterface.hpp)

# Sets install directory for all the headers in the list
foreach(_header ${headers_list})
   install(FILES ${_header} DESTINATION include)
endforeach()
