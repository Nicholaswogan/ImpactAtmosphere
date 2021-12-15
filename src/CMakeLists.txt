
add_subdirectory(minpack)

add_library(impactatmosphere
EvolveAtmFort.f90
HCN_transport.f90
hydrolysis_rates.f90
climate.f90
climate_wrapper.f90)
target_link_libraries(impactatmosphere minpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
