# - Try to find Eigen2 lib
# Once done this will define
#
#  EIGEN_FOUND - system has eigen lib
#  EIGEN_INCLUDE_DIR - the eigen include directory

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

if (EIGEN_INCLUDE_DIRS)
	set(EIGEN_FIND_QUIETLY TRUE)
endif()

find_path(EIGEN_INCLUDE_DIRS NAMES Eigen/Core
	PATHS $ENV{EIGEN_ROOT} ${EIGEN_ROOT}
	PATH_SUFFIXES eigen2
	HINTS
	${INCLUDE_INSTALL_DIR}
	${KDE4_INCLUDE_DIR}
	)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen2 DEFAULT_MSG EIGEN_INCLUDE_DIRS)

mark_as_advanced(EIGEN_INCLUDE_DIRS)


