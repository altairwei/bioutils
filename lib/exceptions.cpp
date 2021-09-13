#include "exceptions.h"

#include "utils.h"

BIOUTILS_BEGIN_SUB_NAMESPACE(utils)

UnknownNucleotideError::UnknownNucleotideError(const char base)
    : std::runtime_error("Unknown nucleotide: '" + characterPrintable(base) + "'")
{

}

std::ostream& operator<<(std::ostream &strm, const UnknownNucleotideError &e) {
  return strm << e.what();
}

BIOUTILS_END_SUB_NAMESPACE(utils)