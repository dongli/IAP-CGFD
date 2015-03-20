set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (leapfrog.exe ${sources})
target_link_libraries (leapfrog.exe geomtk)
