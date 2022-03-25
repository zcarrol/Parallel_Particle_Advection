cmake_minimum_required(VERSION 3.10)

PROJECT(final_proj)

SET(VTK_DIR /home/zacree/Desktop/CMake_SRC/vtk_folder/vtk_build)
file(GLOB SOURCES "*.cxx")

add_executable(final_proj ${SOURCES})

find_package(VTK REQUIRED)

SET(CMAKE_CXX_LINK_FLAGS "-IGL")
if(VTK_LIBRARIES)
  target_link_libraries(final_proj ${VTK_LIBRARIES})
else()
  target_link_libraries(final_proj vtkHybrid)
endif()

vtk_module_autoinit(
  TARGETS final_proj 
  MODULES ${VTK_LIBRARIES})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(final_proj OpenMP::OpenMP_CXX)
endif()
