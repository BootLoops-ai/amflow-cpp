# amflow::ibp sources.

target_sources(amflow PRIVATE
    libp_deriv.cpp
    kira_yaml.cpp
    kira_parse.cpp
    kira_run.cpp
    ibp_cache.cpp
    numd_subst.cpp
    reduce.cpp
)

# numd_subst.cpp streams Kira SYSTEM_*.gz files (AMFLOW_DIFFEQ_NUMD
# numeric-d substitution) via zlib.
find_package(ZLIB REQUIRED)
target_link_libraries(amflow PRIVATE ZLIB::ZLIB)
