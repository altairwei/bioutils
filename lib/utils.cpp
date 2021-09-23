#include "utils.h"

#include <iostream>
#include <bitset>

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

void printHashBit(unsigned long long hash)
{
    std::cout << std::bitset<8*sizeof(hash)>(hash) << std::endl;
}

BIOUTILS_END_SUB_NAMESPACE(utils)