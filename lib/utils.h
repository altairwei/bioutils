#ifndef LIB_UTILS_H
#define LIB_UTILS_H

#include <string>

#include "global.h"

BIOUTILS_BEGIN_SUB_NAMESPACE(utils)

std::string characterPrintable(char c) noexcept(true);
void printHashBit(unsigned long long hash);

BIOUTILS_END_SUB_NAMESPACE(utils)



#endif // LIB_UTILS_H