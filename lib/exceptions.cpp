#include "exceptions.h"

#include "utils.h"

BEGIN_NAMESPACE(bioutils)
BEGIN_NAMESPACE(utils)

UnknownBaseError::UnknownBaseError(const char base)
    : std::runtime_error("Unknown base: '" + characterPrintable(base) + "'")
{

}

std::ostream& operator<<(std::ostream &strm, const UnknownBaseError &e) {
  return strm << e.what();
}

END_NAMESPACE(utils)
END_NAMESPACE(bioutils)