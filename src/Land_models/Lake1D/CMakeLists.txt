# 1-D lake model (module_sf_lake) and its WRF-Hydro driver (module_lake1d).
# Built in every configuration: used by the offline Noah-MP driver and by the
# NUOPC cap.
add_library(hydro_lake1d STATIC
	module_sf_lake.F
	module_lake1d.F
)

add_dependencies(hydro_lake1d hydro_routing hydro_orchestrator hydro_data_rec hydro_mpp noahmp_util)

target_link_libraries(hydro_lake1d PUBLIC
	noahmp_util
	hydro_data_rec
	hydro_orchestrator
	hydro_mpp
	${NETCDF_LIBRARIES}
)
