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

#define REPEAT_LIST_0(x)
#define REPEAT_LIST_1(x) x
#define REPEAT_LIST_2(x) x, REPEAT_LIST_1(x)
#define REPEAT_LIST_3(x) x, REPEAT_LIST_2(x)
#define REPEAT_LIST_4(x) x, REPEAT_LIST_3(x)
#define REPEAT_LIST_5(x) x, REPEAT_LIST_4(x)
#define REPEAT_LIST_6(x) x, REPEAT_LIST_5(x)
#define REPEAT_LIST_7(x) x, REPEAT_LIST_6(x)
#define REPEAT_LIST_8(x) x, REPEAT_LIST_7(x)
#define REPEAT_LIST_9(x) x, REPEAT_LIST_8(x)
#define REPEAT_LIST_10(x) x, REPEAT_LIST_9(x)
#define REPEAT_LIST_20(x) REPEAT_LIST_10(x), REPEAT_LIST_10(x)
#define REPEAT_LIST_30(x) REPEAT_LIST_10(x), REPEAT_LIST_20(x)
#define REPEAT_LIST_40(x) REPEAT_LIST_10(x), REPEAT_LIST_30(x)
#define REPEAT_LIST_50(x) REPEAT_LIST_10(x), REPEAT_LIST_40(x)
#define REPEAT_LIST_60(x) REPEAT_LIST_10(x), REPEAT_LIST_50(x)
#define REPEAT_LIST_70(x) REPEAT_LIST_10(x), REPEAT_LIST_60(x)
#define REPEAT_LIST_80(x) REPEAT_LIST_10(x), REPEAT_LIST_70(x)
#define REPEAT_LIST_90(x) REPEAT_LIST_10(x), REPEAT_LIST_80(x)
#define REPEAT_LIST_100(x) REPEAT_LIST_10(x), REPEAT_LIST_90(x)
#define REPEAT_LIST_N(x, N) REPEAT_LIST_##N(x)

BIOUTILS_END_SUB_NAMESPACE(utils)



#endif // LIB_UTILS_H