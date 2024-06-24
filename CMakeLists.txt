cmake_minimum_required(VERSION 3.14)

project(scran_model_gene_variances
    VERSION 1.0.0
    DESCRIPTION "Aggregate expression values across cells"
    LANGUAGES CXX)

add_library(scran_model_gene_variances INTERFACE)
target_compile_features(scran_model_gene_variances INTERFACE cxx_std_17)
target_include_directories(scran_model_gene_variances INTERFACE include/scran)

add_subdirectory(extern)
target_link_libraries(scran_model_gene_variances INTERFACE tatami::tatami ltla::WeightedLowess scran_core_utils)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()
