add_subdirectory("surfex")

add_library(noahmp_phys STATIC
        module_sf_noahmpdrv.F
        module_sf_noahmp_glacier.F
        module_sf_noahmplsm.F
)

target_link_libraries(noahmp_phys noahmp_submodule snowcro crocus_surfex)
add_dependencies(noahmp_phys noahmp_submodule)

# module_sf_noahmpdrv.F reads the namelist-derived physics options straight
# from config_base's noahmp_io_singleton (see namelist_migration_plan.md).
# CMP0079 NEW is required to link to hydro_orchestrator, defined in a
# different directory (src/OrchestratorLayer).
cmake_policy(SET CMP0079 NEW)
target_link_libraries(noahmp_phys hydro_orchestrator)
add_dependencies(noahmp_phys hydro_orchestrator)
