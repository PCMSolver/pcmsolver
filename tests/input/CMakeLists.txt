set(list_of_libraries utils getkw)

# Python embedding might have been disabled.
# Add tests and references accordingly
if(EMBEDDED_PYTHON)
    list(APPEND external_libraries ${PYTHON_LIBRARIES} ${Boost_PYTHON_LIBRARY})
    add_boosttest(input_embeddedPython "${list_of_libraries}" "${external_libraries}")
    add_reference(tsless.inp ${CMAKE_CURRENT_BINARY_DIR})
    add_reference(restart.inp ${CMAKE_CURRENT_BINARY_DIR})
    add_reference(wavelet.inp ${CMAKE_CURRENT_BINARY_DIR})
else()
    add_boosttest(input "${list_of_libraries}")
    add_reference(@tsless.inp ${CMAKE_CURRENT_BINARY_DIR})
    add_reference(@restart.inp ${CMAKE_CURRENT_BINARY_DIR})
    add_reference(@wavelet.inp ${CMAKE_CURRENT_BINARY_DIR})
endif()
