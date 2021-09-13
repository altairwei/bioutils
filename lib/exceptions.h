#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>
#include <ostream>

#include "global.h"

BIOUTILS_BEGIN_SUB_NAMESPACE(utils)

class UnknownNucleotideError : public std::runtime_error {

public:
    UnknownNucleotideError(const char base);

};

std::ostream& operator<<(std::ostream &strm, const UnknownNucleotideError &e);

BIOUTILS_END_SUB_NAMESPACE(utils)

#endif // EXCEPTIONS_H