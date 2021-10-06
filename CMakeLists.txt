CMAKE_MINIMUM_REQUIRED (VERSION 3.0) 
PROJECT (aneurysm_model) 

# Extends cmake module path - so that FindwaLBerla.cmake in the current directory is found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${aneurysm_model_SOURCE_DIR})
find_package(waLBerla)

add_subdirectory(src)
