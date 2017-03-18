add_executable (ffsl_1d.exe "${CMAKE_CURRENT_SOURCE_DIR}/main_1d.cpp")
target_link_libraries (ffsl_1d.exe geomtk)

add_executable (ffsl_2d.exe "${CMAKE_CURRENT_SOURCE_DIR}/main_2d.cpp")
target_link_libraries (ffsl_2d.exe geomtk)
