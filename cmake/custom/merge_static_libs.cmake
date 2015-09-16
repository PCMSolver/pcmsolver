include(MergeStaticLibs)

get_property(_libraries GLOBAL PROPERTY PCMSolver_LIBRARIES)
merge_static_libs(pcm "${_libraries}")
install(TARGETS pcm ARCHIVE DESTINATION lib)
