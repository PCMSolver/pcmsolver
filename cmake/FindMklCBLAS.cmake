include(CBlasFunctions)
init_vendor_cblas(MKL)
find_cblas_include_dirs(MKL mkl_cblas.h)

if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
find_cblas_libraries(MKL lib/em64t
	mkl_core mkl_intel_lp64 mkl_sequential guide pthread m)
else()
find_cblas_libraries(MKL lib/32 
	mkl_core mkl_intel mkl_sequential guide pthread m)
endif()

if(MKL_CBLAS_LIBRARIES)
	set(MKL_CBLAS_LIBRARIES -Wl,--start-group
		${MKL_CBLAS_LIBRARIES} 
		-Wl,--end-group )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MklCBLAS DEFAULT_MSG
	MKL_CBLAS_INCLUDE_DIRS MKL_CBLAS_LIBRARIES)

mark_as_advanced(MKL_CBLAS_INCLUDE_DIRS MKL_CBLAS_LIBRARIES)
