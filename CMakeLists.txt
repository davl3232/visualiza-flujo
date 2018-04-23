
cmake_minimum_required(VERSION 2.8)
 
PROJECT(ReadOBJ)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(ReadOBJ ReadOBJ.cxx )
 
target_link_libraries(ReadOBJ ${VTK_LIBRARIES})
