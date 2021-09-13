#ifndef DATAIO_H
#define DATAIO_H

#include <string>

#include "global.h"

BIOUTILS_BEGIN_SUB_NAMESPACE(IO)

std::string read_file(const std::string &fileName);
char *read_file(const char *);
char *read_stdin();
std::string read_input(const std::string &argument);

BIOUTILS_END_SUB_NAMESPACE(IO)

#endif //DATAIO_H