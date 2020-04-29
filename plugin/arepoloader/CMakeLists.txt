cmake_minimum_required (VERSION 3.1)
project (arepolibs)
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
find_package(Python3 REQUIRED)
#find_package(OpenVDB REQUIRED)

message(STATUS "PNG package was found: '${PNG_INCLUDE_DIRS}','${PNG_LIBRARIES}'")
message(STATUS "OPENMP package was found:'${OpenMP_CXX_LIBRARIES}'")
message(STATUS "MPI package was found: '${MPI_CXX_INCLUDE_DIRS}','${MPI_CXX_LIBRARIES}'")
message(STATUS "HDF5 package was found: '${HDF5_INCLUDE_DIRS}','${HDF5_LIBRARIES}'")
message(STATUS "GSL package was found: '${GSL_INCLUDE_DIRS}','${GSL_LIBRARIES}'")
message(STATUS "GMP package was found: '${GMP_INCLUDES}','${GMP_LIBRARIES}'")
message(STATUS "PYTHON3 package was found: '${Python3_INCLUDE_DIRS}','${Python3_LIBRARIES}'")
#message(STATUS "OpenVDB package was found: '${OpenVDB_INCLUDE_DIRS}','${OpenVDB_LIBRARIES}'")
#add_executable(illustrisloader 
#main.cpp
#ovdb.cpp

#)
add_library(arepolibs INTERFACE)

target_include_directories(arepolibs INTERFACE 
"${PNG_INCLUDE_DIRS}"
"${MPI_CXX_INCLUDE_DIRS}"
"${HDF5_INCLUDE_DIRS}"
"${GSL_INCLUDE_DIRS}"
"${GMP_INCLUDES}"
"${Python3_INCLUDE_DIRS}"
#"${OpenVDB_INCLUDE_DIRS}"
"${PROJECT_SOURCE_DIR}"
"${PROJECT_SOURCE_DIR}/../../../ArepoVTK/arepo/include/"
"${PROJECT_SOURCE_DIR}/../../../ArepoVTK/include/"

)

target_link_libraries(arepolibs INTERFACE
"${PROJECT_SOURCE_DIR}/../../../ArepoVTK/libarepovtk.so"
"${PROJECT_SOURCE_DIR}/../../../ArepoVTK/arepo/libarepo.so"
"${PNG_LIBRARIES}"
"${OpenMP_CXX_LIBRARIES}"
"${MPI_CXX_LIBRARIES}"
"${HDF5_LIBRARIES}"
"${GSL_LIBRARIES}"
"${GMP_LIBRARIES}"
"${Python3_LIBRARIES}"
pybind11::module
#"${OpenVDB_LIBRARIES}"

#OpenVDB::openvdb

)

find_package(pybind11 REQUIRED)
add_library(arepopy MODULE "py_arepoloader.cpp")
target_link_libraries(arepopy PRIVATE liblm pybind11::module)
set_target_properties(arepopy PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}")
set_target_properties(arepopy PROPERTIES FOLDER "lm/plugin")
set_target_properties(arepopy PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")



include(LmAddPlugin)




lm_add_plugin(
        NAME arepoloader
        LIBRARIES arepolibs
        SOURCES
            "arepoloader.cpp")










