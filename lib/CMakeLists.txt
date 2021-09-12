set(LIBBIOUTILS_HEADERS
    dataio.h
    pattern.h
    exceptions.h
    utils.h
)

set(LIBBIOUTILS_SOURCES
    dataio.cpp
    pattern.cpp
    exceptions.cpp
    utils.cpp
)

add_library(bioutils ${LIBBIOUTILS_HEADERS} ${LIBBIOUTILS_SOURCES})