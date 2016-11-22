if(CMAKE_SYSTEM_NAME MATCHES Windows)
  # Build the library with the __declspec(dllexport) marker
  add_definitions(-DPCMSOLVER_BUILD_SHARED)
  # Get math constants such as M_PI
  add_definitions(-D_USE_MATH_DEFINES)
  # Mark as undefined symbols that have to be defined by the host
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-U,host_writer")
endif()
