#ifndef DATAIO_H
#define DATAIO_H

#include <string>

std::string read_file(const std::string &fileName);
char *read_file(const char *);
char *read_stdin();

#endif //DATAIO_H