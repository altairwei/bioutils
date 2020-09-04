#ifndef LIB_PATTERN_H
#define LIB_PATTERN_H

#include <functional>
#include <string>
#include <set>
#include <vector>

namespace bioutils {

namespace algorithms {

#define ISNTP(C) ((C) == 'A' || (C) == 'T' || (C) == 'C' || (C) == 'G' \
        || (C) == 'a' || (C) == 't' || (C) == 'c' || (C) == 'g')

typedef unsigned long long hash_t;

hash_t hash_kmer(const char *, size_t);
bool is_ntp(char c);
int bpton(char);

enum PatternCountAlgorithms { BruteForce, BruteForceByHand, RabinKarp };

void find_do(const char *text, const char *pattern,
    std::function<void(const size_t, const char *, const char *)> callback);
std::size_t PatternCount(const char *text, const char *pattern, enum PatternCountAlgorithms algo);
void FrequentWords(const std::string text, const int k, std::set<std::string> &output);

} // algorithms

} // bioutils

#endif // LIB_PATTERN_H