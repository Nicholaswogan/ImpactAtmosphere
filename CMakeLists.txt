cmake_minimum_required(VERSION "3.12")
project(impactatmosphere Fortran C)

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

option(CMAKE_POSITION_INDEPENDENT_CODE "Makes code position independent." ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" 
   AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wunused")
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_subdirectory(src)

add_executable(test src/test.f90)
target_link_libraries(test impactatmosphere)

# python extension
set(F2PY_SOURCES 
${CMAKE_SOURCE_DIR}/src/EvolveAtmFort.f90 
${CMAKE_SOURCE_DIR}/src/HCN_transport.f90 
${CMAKE_SOURCE_DIR}/src/hydrolysis_rates.f90
${CMAKE_SOURCE_DIR}/src/climate_wrapper.f90)

if (NOT SKBUILD)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
endif()
  
find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)
find_package(F2PY REQUIRED)
  
message(STATUS "F2PY included from: ${F2PY_INCLUDE_DIR}")
  
set(f2py_module_name "EvolveAtmFort")
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
  
  add_custom_target(${f2py_module_name} ALL
  DEPENDS ${F2PY_SOURCES}
  )
  
add_custom_command(
OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
COMMAND ${F2PY_EXECUTABLE}
  -m ${f2py_module_name}
  ${F2PY_SOURCES}
  --lower
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
DEPENDS ${F2PY_SOURCES}
)
  
add_library(${generated_module_file} MODULE
            "${f2py_module_name}module.c"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${f2py_module_name}-f2pywrappers2.f90"
             ${PHOTOCHEM_SOURCES})
    
target_link_libraries(${generated_module_file} impactatmosphere)
target_include_directories(${generated_module_file} PUBLIC
                          "${CMAKE_BINARY_DIR}/modules"
                          ${F2PY_INCLUDE_DIRS}
                          ${PYTHON_INCLUDE_DIRS})                            
set_target_properties(${generated_module_file} PROPERTIES SUFFIX "")
set_target_properties(${generated_module_file} PROPERTIES PREFIX "")    
  
if (UNIX)
  if (APPLE)
    set_target_properties(${generated_module_file} PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${generated_module_file} PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()
  
if (SKBUILD)
  install(TARGETS ${generated_module_file} DESTINATION ImpactAtmosphere)
else()
  install(TARGETS ${generated_module_file} DESTINATION ${CMAKE_SOURCE_DIR}/ImpactAtmosphere)
endif()
