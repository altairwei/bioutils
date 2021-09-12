#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>
#include <ostream>

#include "global.h"

BEGIN_NAMESPACE(bioutils)
BEGIN_NAMESPACE(utils)

class UnknownBaseError : public std::runtime_error {

public:
    UnknownBaseError(const char base);

};

std::ostream& operator<<(std::ostream &strm, const UnknownBaseError &e);

END_NAMESPACE(utils)
END_NAMESPACE(bioutils)

#endif // EXCEPTIONS_H