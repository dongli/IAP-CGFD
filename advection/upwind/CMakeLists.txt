set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (upwind.exe ${sources})
target_link_libraries (upwind.exe geomtk)
