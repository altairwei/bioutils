cmake_minimum_required(VERSION 3.0)

project(bioutils)

option(BUILD_TESTING "Building the test target or not." ON)

set(EXE_PREFIX "bio")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CONAN_INSTALL_MANUALLY)
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
        message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
        file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/master/conan.cmake"
                        "${CMAKE_BINARY_DIR}/conan.cmake")
    endif()
    include(${CMAKE_BINARY_DIR}/conan.cmake)
    conan_check(REQUIRED)
    conan_cmake_run(
        CONANFILE conanfile.py
        INSTALL_FOLDER ${CMAKE_BINARY_DIR}
        BUILD_TYPE ${CMAKE_BUILD_TYPE}
        BUILD missing
    )
endif()

# set output path
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_output_dirs_setup()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_BINARY_DIR})
find_package(argparse REQUIRED)
find_package(gtest REQUIRED)

add_subdirectory(lib)
add_subdirectory(src/find)
add_subdirectory(src/seq)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()