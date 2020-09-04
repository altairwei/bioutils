set(LIBBIOUTILS_HEADERS
    dataio.h
    pattern.h
)

set(LIBBIOUTILS_SOURCES
    dataio.cpp
    pattern.cpp
)

add_library(bioutils ${LIBBIOUTILS_HEADERS} ${LIBBIOUTILS_SOURCES})