if(CMAKE_SYSTEM_NAME MATCHES Windows)
  # Get math constants such as M_PI
  add_definitions(-D_USE_MATH_DEFINES)
endif()
