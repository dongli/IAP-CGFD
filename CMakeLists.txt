cmake_minimum_required (VERSION 2.8)

project (IAP-CGFD CXX)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# GEOMTK
add_subdirectory ("${PROJECT_SOURCE_DIR}/external/geomtk")
include_directories (${GEOMTK_INCLUDE_DIRS})
link_directories (${GEOMTK_LIBRARY_DIRS})

add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/upwind")
