# Per-domain source files.
#
# Each domain populates target_sources(amflow PRIVATE ...).
# Order follows the dependency DAG documented in docs/ARCHITECTURE.md:

add_subdirectory(numeric)
add_subdirectory(algebra)
add_subdirectory(ode)
add_subdirectory(qft)
add_subdirectory(ibp)
add_subdirectory(pipeline)
add_subdirectory(api)
add_subdirectory(cli)
