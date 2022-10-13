cmake_minimum_required(VERSION "3.14")
project(impactatmosphere LANGUAGES Fortran C)
include(cmake/CPM.cmake)

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

option(CMAKE_POSITION_INDEPENDENT_CODE "Makes code position independent." ON)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" 
   AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wunused")
endif()

if (NOT SKBUILD)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
endif()

find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)
find_package(F2PY REQUIRED)
find_package(Cython REQUIRED)

set(f2py_module_name "EvolveAtmFort")
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
add_subdirectory(src)

# test
add_executable(test src/test.f90)
target_link_libraries(test impactatmosphere)
  
# install
if (SKBUILD)
  install(TARGETS ${generated_module_file} _cvode DESTINATION ImpactAtmosphere)
else()
  install(TARGETS ${generated_module_file} _cvode DESTINATION ${CMAKE_SOURCE_DIR}/ImpactAtmosphere)
endif()