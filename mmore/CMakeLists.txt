cmake_minimum_required(VERSION 3.10)
include(ExternalProject)

project(mmoreseqs VERSION 1.0.0)

set(CMAKE_C_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY build)

# -----------
# Executables
# -----------

add_executable(
        mmoreseqs
        src/application.c
)

# ------------
# Dependencies
# ------------

ExternalProject_Add(
        easel
        GIT_REPOSITORY git@github.com:EddyRivasLab/easel.git
        GIT_TAG Bio-Easel-0.14
        PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib/easel
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoconf
        COMMAND ./configure --prefix=${SOURCE_DIR}
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        TEST_COMMAND ""
)
ExternalProject_Get_property(easel SOURCE_DIR)
target_include_directories(mmoreseqs PRIVATE ${SOURCE_DIR})
add_dependencies(mmoreseqs easel)
