add_library(hydro_routing_rivertemperature STATIC
	module_channel_routing_local.F
	module_lsm_forcing_local.F
	module_river_temperature.F
	module_sf_noahmpdrv_local.F
        weather_disaggrgate.F
)

add_dependencies(hydro_routing_rivertemperature hydro_utils)
