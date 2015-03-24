set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (beam_warming.exe ${sources})
target_link_libraries (beam_warming.exe geomtk)
