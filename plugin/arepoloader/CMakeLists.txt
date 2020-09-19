cmake_minimum_required (VERSION 3.1)
project (arepoloader)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_C_FLAGS -DH5_USE_110_API)

set( CMAKE_BUILD_TYPE Release)
add_compile_options(-O2)
#add_compile_options(-g)
#placed findGMP.cmake in this directory...
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} )#${PROJECT_SOURCE_DIR}/../../openvdb/cmake/)

find_package(PNG REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)
find_package(GSL REQUIRED)
find_package(GMP REQUIRED)

message(STATUS "PNG package was found: '${PNG_INCLUDE_DIRS}','${PNG_LIBRARIES}'")
message(STATUS "OPENMP package was found:'${OpenMP_CXX_LIBRARIES}'")
message(STATUS "MPI package was found: '${MPI_CXX_INCLUDE_DIRS}','${MPI_CXX_LIBRARIES}'")
message(STATUS "HDF5 package was found: '${HDF5_INCLUDE_DIRS}','${HDF5_LIBRARIES}'")
message(STATUS "GSL package was found: '${GSL_INCLUDE_DIRS}','${GSL_LIBRARIES}'")
message(STATUS "GMP package was found: '${GMP_INCLUDES}','${GMP_LIBRARIES}'")


add_library(arepoloader SHARED "pluecker.h" "statstags.h" "arepoloader.h" "raytetrapluecker.cpp" "arepoloader.cpp" "mesh_arepo.cpp" "medium_het_arepo.cpp" "renderer_volpt_arepo.cpp" "renderer_arepo_vrls.cpp" "lightscene.cpp")

option(MOCK_AREPO "Mock Arepo Data" OFF)
if(MOCK_AREPO)
    message("use arepo data mockup (only one tetrahedron)")
    target_compile_definitions(arepoloader PUBLIC MOCK_AREPO)
endif()


target_include_directories(arepoloader PUBLIC 
"${PNG_INCLUDE_DIRS}"
"${MPI_CXX_INCLUDE_DIRS}"
"${HDF5_INCLUDE_DIRS}"
"${GSL_INCLUDE_DIRS}"
"${GMP_INCLUDES}"
"/home/beisel/ArepoVTK/arepo/include/"
"/home/beisel/ArepoVTK/include/"
)

target_link_libraries(arepoloader PUBLIC
"/home/beisel/ArepoVTK/libarepovtk.a"
"/home/beisel/ArepoVTK/arepo/libarepo.a"
"${PROJECT_SOURCE_DIR}"
"${PNG_LIBRARIES}"
"${OpenMP_CXX_LIBRARIES}"
"${MPI_CXX_LIBRARIES}"
"${HDF5_LIBRARIES}"
"${GSL_LIBRARIES}"
"${GMP_LIBRARIES}"
"${Python3_LIBRARIES}"
pybind11::module
liblm
)


target_compile_definitions(arepoloader PRIVATE -DLM_AREPO_EXPORTS)
set_target_properties(arepoloader PROPERTIES PREFIX "")
set_target_properties(arepoloader PROPERTIES DEBUG_POSTFIX "-debug")
set_target_properties(arepoloader PROPERTIES FOLDER "lm/plugin")
set_target_properties(arepoloader PROPERTIES
ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")




find_package(pybind11 REQUIRED)
add_library(arepopy MODULE "py_arepoloader.cpp")
target_link_libraries(arepopy PRIVATE liblm pybind11::module arepoloader)
set_target_properties(arepopy PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}")
set_target_properties(arepopy PROPERTIES FOLDER "lm/plugin")
set_target_properties(arepopy PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


#dont needed?
#include(LmAddPlugin)
#lm_add_plugin(
#        NAME arepoloader
#        LIBRARIES arepolibs
#        SOURCES
#            "arepoloader.cpp")










