set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (semi_lagrangian.exe ${sources})
target_link_libraries (semi_lagrangian.exe geomtk)
