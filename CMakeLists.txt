cmake_minimum_required (VERSION 2.8)

project (IAP-CGFD CXX)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DBOOST_VARIANT_USE_RELAXED_GET_BY_DEFAULT")

# GEOMTK
add_subdirectory ("${PROJECT_SOURCE_DIR}/external/geomtk")
include_directories (${GEOMTK_INCLUDE_DIRS})
link_directories (${GEOMTK_LIBRARY_DIRS})

# GSL
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
find_package (GSL REQUIRED)
include_directories (${GSL_INCLUDE_DIRS})
link_directories (${GSL_LIBRARY_DIRS})

add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/upwind")
add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/beam_warming")
add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/lax_wendroff")
add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/fromm")
add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/ftcs")
add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/leapfrog")
add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/tspas")
add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/mpdata")
add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/crank_nicolson")
add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/semi_lagrangian")
add_subdirectory ("${PROJECT_SOURCE_DIR}/advection/ffsl")

add_subdirectory ("${PROJECT_SOURCE_DIR}/shallow_water/cartesian/semi_lagrangian")
