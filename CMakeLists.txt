cmake_minimum_required(VERSION 3.0)

project(bioutils)

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
conan_basic_setup()

set(BPFIND_HEADERS
    src/dataio.h
    src/bpfind.h
)

set(BPFIND_SOURCES
    src/dataio.cpp
    src/bpfind.cpp
)

set(BPSEQ_HEADERS
    src/dataio.h
)

set(BPSEQ_SOURCES
    src/dataio.cpp
    src/bpseq.cpp
)

add_executable(bpfind ${BPFIND_HEADERS} ${BPFIND_SOURCES})
add_executable(bpseq ${BPSEQ_HEADERS} ${BPSEQ_SOURCES})
target_link_libraries(bpfind ${CONAN_LIBS})
target_link_libraries(bpseq ${CONAN_LIBS})