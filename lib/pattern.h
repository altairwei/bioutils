#ifndef LIB_PATTERN_H
#define LIB_PATTERN_H

#include <functional>
#include <string>
#include <set>
#include <map>
#include <vector>

#include "global.h"

BIOUTILS_BEGIN_SUB_NAMESPACE(algorithms)

#define ISNTP(C) ((C) == 'A' || (C) == 'T' || (C) == 'C' || (C) == 'G' \
        || (C) == 'a' || (C) == 't' || (C) == 'c' || (C) == 'g')

typedef unsigned long long hash_t;

hash_t PatternToNumber(const std::string &kmer);
bool is_ntp(char c);
int NucleobaseToInt(char);

enum class PatternCountAlgorithms { BruteForce, BruteForceByHand, RabinKarp };

void find_do(const char *text, const char *pattern,
    std::function<void(const size_t, const char *, const char *)> callback);
std::size_t PatternCount(const char *text, const char *pattern, enum PatternCountAlgorithms algo);
std::vector<size_t> PatternIndex(const char *text, const char *pattern);

enum class AlgorithmEfficiency {Slow, Fast, Faster, Fastest};
std::set<std::string> FrequentWords(const std::string &text, const int k, AlgorithmEfficiency algo = AlgorithmEfficiency::Slow);
std::set<std::string> FrequentWordsSlow(const std::string &text, const int k);
std::set<std::string> FrequentWordsFast(const std::string &text, const int k);
std::map<std::string, size_t> FrequencyTable(const std::string &text, const int k);
size_t MaxMap(const std::map<std::string, size_t> &input_map);

BIOUTILS_END_SUB_NAMESPACE(algorithms)

#endif // LIB_PATTERN_H