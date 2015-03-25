set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (tspas.exe ${sources})
target_link_libraries (tspas.exe geomtk)
