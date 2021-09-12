#include "exceptions.h"

#include "utils.h"

BEGIN_NAMESPACE(bioutils)
BEGIN_NAMESPACE(utils)

UnknownNucleotideError::UnknownNucleotideError(const char base)
    : std::runtime_error("Unknown nucleotide: '" + characterPrintable(base) + "'")
{

}

std::ostream& operator<<(std::ostream &strm, const UnknownNucleotideError &e) {
  return strm << e.what();
}

END_NAMESPACE(utils)
END_NAMESPACE(bioutils)