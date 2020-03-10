cmake_minimum_required(VERSION 3.5)
project(gridspec VERSION 0.0.0 LANGUAGES C)
 
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMake)
find_package(NetCDF REQUIRED)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src)