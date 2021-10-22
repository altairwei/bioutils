#ifndef LIB_UTILS_H
#define LIB_UTILS_H

#include <string>
#include <stdexcept>
#include <vector>

#include "global.h"

BIOUTILS_BEGIN_SUB_NAMESPACE(utils)

std::string characterPrintable(char c) noexcept(true);
void printHashBit(unsigned long long hash);

/*!
    Find the maximum value of a \a input_map
 */
template <typename T>
size_t MaxMap(const T &input_map) noexcept(false)
{
    if (input_map.empty())
        throw std::runtime_error("input map is empty.");

    auto p = input_map.cbegin();
    size_t max = (*p++).second;

    while (p != input_map.cend()) {
        size_t cur = (*p++).second;
        if (cur > max) max = cur;
    }

    return max;
}

template <typename T>
size_t MaxArray(const T &input_array) noexcept(false)
{
    if (input_array.empty())
        throw std::runtime_error("input array is empty.");

    auto p = input_array.cbegin();
    size_t max = (*p++);

    while (p != input_array.cend()) {
        size_t cur = (*p++);
        if (cur > max) max = cur;
    }

    return max;
}

BIOUTILS_END_SUB_NAMESPACE(utils)



#endif // LIB_UTILS_H