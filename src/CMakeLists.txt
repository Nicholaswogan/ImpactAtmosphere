
add_subdirectory(minpack)
add_subdirectory(lapack)

add_library(impactatmosphere
EvolveAtmFort.f90
HCN_transport.f90
hydrolysis_rates.f90
climate.f90
climate_wrapper.f90)
target_link_libraries(impactatmosphere minpack part_of_lapack)

# python extension
set(F2PY_SOURCES 
${CMAKE_CURRENT_SOURCE_DIR}/EvolveAtmFort.f90 
${CMAKE_CURRENT_SOURCE_DIR}/HCN_transport.f90 
${CMAKE_CURRENT_SOURCE_DIR}/hydrolysis_rates.f90
${CMAKE_CURRENT_SOURCE_DIR}/climate_wrapper.f90)
  
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
            "${f2py_module_name}-f2pywrappers2.f90")
    
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

add_subdirectory(cvode)
