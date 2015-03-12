set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (lax_wendroff.exe ${sources})
target_link_libraries (lax_wendroff.exe geomtk)
