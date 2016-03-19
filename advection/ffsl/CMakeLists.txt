set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (ffsl.exe ${sources})
target_link_libraries (ffsl.exe geomtk)
