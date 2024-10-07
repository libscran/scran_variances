cmake_minimum_required(VERSION 3.14)

project(scran_variances
    VERSION 0.1.2
    DESCRIPTION "Model per-gene variances in single-cell expression"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Library
add_library(scran_variances INTERFACE)
add_library(libscran::scran_variances ALIAS scran_variances)

target_include_directories(scran_variances INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/scran_variances>)
target_compile_features(scran_variances INTERFACE cxx_std_17)

# Dependencies
option(SCRAN_VARIANCES_FETCH_EXTERN "Automatically fetch scran_variances's external dependencies." ON)
if(SCRAN_VARIANCES_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(tatami_tatami 3.0.0 CONFIG REQUIRED)
    find_package(tatami_tatami_stats 1.0.0 CONFIG REQUIRED)
    find_package(ltla_WeightedLowess 2.0.0 CONFIG REQUIRED)
    find_package(libscran_scran_blocks 0.1.0 CONFIG REQUIRED)
endif()

target_link_libraries(scran_variances INTERFACE tatami::tatami tatami::tatami_stats ltla::WeightedLowess libscran::scran_blocks)

# Tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(SCRAN_VARIANCES_TESTS "Build scran_variances's test suite." ON)
else()
    option(SCRAN_VARIANCES_TESTS "Build scran_variances's test suite." OFF)
endif()

if(SCRAN_VARIANCES_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif() 
endif()

# Install
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/scran_variances)

install(TARGETS scran_variances
    EXPORT scran_variancesTargets)

install(EXPORT scran_variancesTargets
    FILE libscran_scran_variancesTargets.cmake
    NAMESPACE libscran::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_scran_variances)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_variancesConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_scran_variances)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_variancesConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_variancesConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_variancesConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_scran_variances)
