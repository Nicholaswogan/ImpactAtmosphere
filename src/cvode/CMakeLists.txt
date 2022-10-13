CPMAddPackage(
  NAME cvode 
  VERSION 5.7.0
  OPTIONS 
    "BUILD_SHARED_LIBS OFF" 
    "SUNDIALS_INDEX_SIZE 32"
    "SUNDIALS_PRECISION double"
  URL "https://github.com/LLNL/sundials/releases/download/v5.7.0/cvode-5.7.0.tar.gz"
  EXCLUDE_FROM_ALL ON
)

add_cython_target(_cvode _cvode.pyx)
add_library(_cvode MODULE ${_cvode})
python_extension_module(_cvode)
target_link_libraries(_cvode sundials_cvode)
target_include_directories(_cvode PUBLIC ${NumPy_INCLUDE_DIRS})

