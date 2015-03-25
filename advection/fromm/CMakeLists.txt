set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (fromm.exe ${sources})
target_link_libraries (fromm.exe geomtk)
