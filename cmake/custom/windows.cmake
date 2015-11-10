set(EXE "")
if(CMAKE_SYSTEM_NAME MATCHES Windows)
  # Build the library with the __declspec(dllexport) marker
  add_definitions(-DPCMSOLVER_BUILD_SHARED)
  # Get math constants such as M_PI
  add_definitions(-D_USE_MATH_DEFINES)
  set(EXE ".exe")
endif()
