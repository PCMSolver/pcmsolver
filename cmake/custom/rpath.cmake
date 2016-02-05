if(APPLE)
    if(NOT DEFINED CMAKE_MACOSX_RPATH)
        set(CMAKE_MACOSX_RPATH ON)
        # Mark as undefined symbols that have to be defined by the host
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-U,_host_writer")
    endif()
endif()
