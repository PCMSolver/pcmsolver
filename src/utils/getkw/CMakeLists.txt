target_sources(pcm-objlib
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Getkw.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GetkwError.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Section.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/messages.cpp
  )

# List of headers
list(APPEND _headers_list
  Getkw.h
  Section.h
  Keyword.h
  GetkwError.h
  )
# Sets install directory for all the headers in the list
foreach(_header IN LISTS _headers_list)
   install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/utils/getkw)
endforeach()
