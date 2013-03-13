include(CBlasFunctions)

init_vendor_cblas(ATLAS)
find_cblas_include_dirs(ATLAS cblas.h)
find_cblas_libraries(ATLAS lib cblas atlas f77blas)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AtlasCBLAS DEFAULT_MSG 
	ATLAS_CBLAS_INCLUDE_DIRS ATLAS_CBLAS_LIBRARIES)

mark_as_advanced(ATLAS_CBLAS_INCLUDE_DIRS ATLAS_CBLAS_LIBRARIES)

