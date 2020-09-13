#ifndef DATAIO_H
#define DATAIO_H

#include <string>

namespace bioutils {

namespace IO {

std::string read_file(const std::string &fileName);
char *read_file(const char *);
char *read_stdin();
std::string read_input(const std::string &argument);

} // IO

} // bioutils

#endif //DATAIO_H