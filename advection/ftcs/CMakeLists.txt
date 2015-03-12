set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (ftcs.exe ${sources})
target_link_libraries (ftcs.exe geomtk)
