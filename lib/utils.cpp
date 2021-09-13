#include "utils.h"

BIOUTILS_BEGIN_SUB_NAMESPACE(utils)

std::string characterPrintable(char c) noexcept(true)
{
    std::string repr;

    switch (c)
    {
    case '\n':
        repr = "\\n";
        break;
    case '\r':
        repr = "\\r";
    case '\f':
        repr = "\\f";
    case '\t':
        repr = "\\t";
    case ' ':
        repr = "whitespace";
    case '\v':
        repr = "\\v";
    default:
        repr = c;
        break;
    }

    return repr;
}

BIOUTILS_END_SUB_NAMESPACE(utils)