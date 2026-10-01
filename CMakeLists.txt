cmake_minimum_required(VERSION 3.16)
project(AMFlowCpp
    VERSION 1.1.0
    DESCRIPTION "AMFlow.cpp — C++17 reimplementation of the AMFlow algorithm (numeric/algebra/ode/qft/ibp/pipeline/api/cli)"
    HOMEPAGE_URL "https://github.com/chang18/amflow-cpp"
    LANGUAGES CXX
)

# ---------------------------------------------------------------------------
# Standards & defaults
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(AMFLOW_BUILD_DRIVER   "Build standalone JSON driver"  ON)
option(AMFLOW_BUILD_EXAMPLES "Build small example programs"  OFF)
option(AMFLOW_INSTALL        "Generate install / export rules" ON)
option(AMFLOW_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
find_package(PkgConfig REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint>=3.0)
# Some symbols we use (e.g. fmpz_mpoly_q_set_str_pretty,
# fmpz_mpoly_q_evaluate_acb) require FLINT 3.x.  Ubuntu 22.04 ships
# FLINT 2.8.4 in libflint-dev — that won't work; build FLINT from
# source or use a distro / Conda channel with FLINT >= 3.0.

if(NOT TARGET nlohmann_json::nlohmann_json)
    find_package(nlohmann_json 3.7.0 QUIET)
    if(NOT nlohmann_json_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            nlohmann_json
            GIT_REPOSITORY https://github.com/nlohmann/json.git
            GIT_TAG        v3.11.3
            GIT_SHALLOW    TRUE
        )
        set(JSON_BuildTests OFF CACHE INTERNAL "")
        # amflow links nlohmann_json PUBLIC (the api/ headers expose
        # nlohmann::json), so the fetched target ends up in the
        # AMFlowCppTargets install export.  Installing the fetched copy
        # alongside puts it in its own export set — otherwise the generate
        # step fails with "requires target 'nlohmann_json' that is not in
        # any export set".  The installed AMFlowCppConfig.cmake then finds
        # it via find_dependency(nlohmann_json).
        set(JSON_Install    ON  CACHE INTERNAL "")
        FetchContent_MakeAvailable(nlohmann_json)
    endif()
endif()

# ---------------------------------------------------------------------------
# Library: amflow  (domain-oriented; sources added under src/<domain>/)
# ---------------------------------------------------------------------------
add_library(amflow STATIC
    src/architecture_marker.cpp
)

target_include_directories(amflow
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

find_package(Threads REQUIRED)

target_link_libraries(amflow
    PUBLIC
        PkgConfig::FLINT
        nlohmann_json::nlohmann_json
        Threads::Threads
)

target_compile_features(amflow PUBLIC cxx_std_17)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(amflow PRIVATE
        -Wall -Wextra -Wpedantic
        -Wno-unused-parameter
    )
    if(AMFLOW_WARNINGS_AS_ERRORS)
        target_compile_options(amflow PRIVATE -Werror)
    endif()
endif()

# Per-domain source files + cli executable.
add_subdirectory(src)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
include(CTest)  # provides the standard BUILD_TESTING option (default ON)
if(BUILD_TESTING)
    pkg_check_modules(GTEST REQUIRED IMPORTED_TARGET gtest)
    add_subdirectory(tests)
endif()

# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------
if(AMFLOW_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------
message(STATUS "")
message(STATUS "AMFlow.cpp ${PROJECT_VERSION}")
message(STATUS "  Build type      : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard    : ${CMAKE_CXX_STANDARD}")
message(STATUS "  FLINT version   : ${FLINT_VERSION}")
message(STATUS "  Tests           : ${BUILD_TESTING}")
message(STATUS "  Driver          : ${AMFLOW_BUILD_DRIVER}")
message(STATUS "  Examples        : ${AMFLOW_BUILD_EXAMPLES}")
message(STATUS "  Install rules   : ${AMFLOW_INSTALL}")
message(STATUS "")

# ---------------------------------------------------------------------------
# Install / export rules
# ---------------------------------------------------------------------------
if(AMFLOW_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    set(AMFLOWCPP_CMAKE_CONFIG_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/AMFlowCpp")

    install(TARGETS amflow
        EXPORT AMFlowCppTargets
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    install(DIRECTORY include/amflow
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp"
    )

    install(EXPORT AMFlowCppTargets
        FILE AMFlowCppTargets.cmake
        NAMESPACE AMFlowCpp::
        DESTINATION ${AMFLOWCPP_CMAKE_CONFIG_DIR}
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/AMFlowCppConfigVersion.cmake"
        VERSION       ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/AMFlowCppConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/AMFlowCppConfig.cmake"
        INSTALL_DESTINATION ${AMFLOWCPP_CMAKE_CONFIG_DIR}
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/AMFlowCppConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/AMFlowCppConfigVersion.cmake"
        DESTINATION ${AMFLOWCPP_CMAKE_CONFIG_DIR}
    )
endif()
