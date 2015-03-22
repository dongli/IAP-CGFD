set (sources
     "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

add_executable (crank_nicolson.exe ${sources})
target_link_libraries (crank_nicolson.exe geomtk gsl gslcblas)
