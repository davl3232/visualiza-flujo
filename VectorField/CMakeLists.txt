
cmake_minimum_required(VERSION 2.8)
 
PROJECT(VectorField)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(VectorField VectorField.cxx )
 
target_link_libraries(VectorField ${VTK_LIBRARIES})
