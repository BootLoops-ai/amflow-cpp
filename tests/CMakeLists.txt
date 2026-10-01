include(GoogleTest)

set(AMFLOW_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_numeric_options.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_numeric_acb_wrap.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_numeric_rational.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_algebra_mpoly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_algebra_matrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_algebra_context_migration.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ode_blocks.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ode_sparse.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ode_asy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ode_regular.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ode_inf.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ode_qqbar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ode_jordan.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ode_acb_rational.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ode_normalize.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ode_zero.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_qft_jintegral.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_qft_family_config.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_qft_family_uf.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_qft_topology.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_qft_amfmode.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_qft_region.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_qft_boundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_qft_vacuum.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ibp_libp_deriv.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ibp_kira.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ibp_kira_stall.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ibp_cache.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ibp_reduce.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_amflow_factorize.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_amflow_amfsystem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_amflow_solve_integrals.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_api_run_json.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_cli_integration.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_amflow_amfsystem_mma_ref.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_amflow_solve_integrals_mma_ref.cpp
)

add_executable(amflow_tests ${AMFLOW_TEST_SOURCES})

target_link_libraries(amflow_tests
    PRIVATE
        amflow
        PkgConfig::GTEST
        pthread
)

add_dependencies(amflow_tests amflow_cli)

target_compile_features(amflow_tests PRIVATE cxx_std_17)
target_compile_definitions(amflow_tests PRIVATE
    AMFLOW_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    AMFLOW_CLI_PATH="$<TARGET_FILE:amflow_cli>"
)

gtest_discover_tests(amflow_tests)
