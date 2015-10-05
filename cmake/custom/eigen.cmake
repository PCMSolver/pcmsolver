set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external/eigen3/include/eigen3)
message(STATUS "Eigen 3.2.0 is located here: " ${EIGEN3_INCLUDE_DIR})
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

if(ENABLE_EIGEN_MKL)
   message(STATUS "ENABLE_EIGEN_MKL option requires Intel MKL 10.3")
   message(STATUS "   Be sure you have read http://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html")
   set(EIGEN_USE_MKL_ALL ON)
endif()
