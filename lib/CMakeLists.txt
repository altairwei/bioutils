set(LIBBIOUTILS_HEADERS
    dataio.h
)

set(LIBBIOUTILS_SOURCES
    dataio.cpp
)

add_library(bioutils ${LIBBIOUTILS_HEADERS} ${LIBBIOUTILS_SOURCES})