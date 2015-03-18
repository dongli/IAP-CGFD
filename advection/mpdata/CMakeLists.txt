set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (mpdata.exe ${sources})
target_link_libraries (mpdata.exe geomtk)
