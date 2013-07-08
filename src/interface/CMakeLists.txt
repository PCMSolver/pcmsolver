set(interface_SOURCES 
        Interface.cpp
)

add_library(interface ${interface_SOURCES})
add_dependencies(interface libgetkw)

install(TARGETS interface ARCHIVE DESTINATION lib)

install(FILES Interface.hpp DESTINATION include)
