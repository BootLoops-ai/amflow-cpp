# amflow_cli — JSON driver.

add_executable(amflow_cli
    main.cpp
)

target_link_libraries(amflow_cli
    PRIVATE
        amflow
)

target_compile_features(amflow_cli PRIVATE cxx_std_17)

if(AMFLOW_INSTALL)
    include(GNUInstallDirs)
    # Preserve the linker-derived RUNPATH (e.g. /usr/local/lib so the
    # installed binary still finds libflint / libmpfr / libgmp), instead
    # of letting CMake rewrite it at install time — the rewrite hits a
    # known mismatch when the build RUNPATH already equals the target
    # install path.
    set_target_properties(amflow_cli PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
    install(TARGETS amflow_cli
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
