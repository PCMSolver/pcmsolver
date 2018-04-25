target_sources(pcm-objlib
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Input.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Meddle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PhysicalConstants.cpp
  )

# List of headers
list(APPEND headers_list
  Input.hpp
  Meddle.hpp
  )
# Sets install directory for all the headers in the list
foreach(_header IN LISTS headers_list)
    install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/interface)
endforeach()
