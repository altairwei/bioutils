#ifndef LIB_PATTERN_H
#define LIB_PATTERN_H

#include <functional>
#include <string>
#include <string_view>
#include <set>
#include <map>
#include <vector>

#include "global.h"

BIOUTILS_BEGIN_SUB_NAMESPACE(algorithms)

#define ISNTP(C) ((C) == 'A' || (C) == 'T' || (C) == 'C' || (C) == 'G' \
        || (C) == 'a' || (C) == 't' || (C) == 'c' || (C) == 'g')

typedef unsigned long long hash_t;
extern const int MAX_HASHABLE_LENGTH;

enum class PatternCountAlgorithms { BruteForce, BruteForceByHand, RabinKarp };
enum class AlgorithmEfficiency {Default, Slow, Fast, Faster, Fastest};

hash_t PatternToNumber(const std::string_view pattern, AlgorithmEfficiency algo = AlgorithmEfficiency::Slow);
hash_t PatternToNumberBitwise(const std::string_view pattern);
hash_t PatternToNumberRecursive(const std::string_view pattern);
std::string NumberToPatternBitwise(const hash_t number, const int length);
bool is_ntp(char c);
std::string ReverseComplement(const std::string_view oriSeq);

void find_do(const std::string_view text, const std::string_view pattern,
    std::function<void(const size_t, const std::string_view, const std::string_view)> callback);
size_t PatternCount(const std::string_view text, const std::string_view pattern, AlgorithmEfficiency algo);
size_t PatternCount_RK(const std::string_view text, const std::string_view pattern);
size_t PatternCount_BF(const std::string_view text, const std::string_view pattern);
std::vector<size_t> PatternIndex(const std::string_view text, const std::string_view pattern);
std::vector<size_t> PatternIndexApproximate(const std::string_view text, const std::string_view pattern, const size_t d);

std::set<std::string> FrequentWords(const std::string_view text, const int k, AlgorithmEfficiency algo = AlgorithmEfficiency::Slow);
std::set<std::string> FrequentWordsSlow(const std::string_view text, const int k);
std::set<std::string> FrequentWordsByPerfectHash(const std::string_view text, const int k);
std::set<std::string> FrequentWordsByStdHash(const std::string_view text, const int k);
std::set<std::string> FrequentWordsBySorting(const std::string_view text, const int k);
std::unordered_map<std::string, uint> FrequencyTable(const std::string_view text, const int k);
std::unordered_map<std::string, uint> FrequencyTableWithMismatches(
    const std::string_view text, const int k, const int d, bool rev_comp = false);
std::vector<uint> FrequencyArray(const std::string_view text, const int k);
std::set<std::string> FrequentWordsWithMismatches(
    const std::string_view text, const int k, const int d, bool rev_comp = false);
std::set<std::string> FrequentWordsWithMismatchesBySorting(
    const std::string_view text, const int k, const int d, bool rev_comp = false);

std::set<std::string> FindClumps(const std::string_view genome, int k, int window_length, int times, AlgorithmEfficiency algo = AlgorithmEfficiency::Default);
std::set<std::string> FindClumpsBetterWithStdHash(const std::string_view genome, int k, int window_length, int times);
std::set<std::string> FindClumpsBetterWithPerfectHash(const std::string_view genome, int k, int window_length, int times);
std::vector<size_t> FindMinimumSkew(const std::string_view genome);

size_t HammingDistance(const std::string_view pattern1, const std::string_view pattern2);
std::set<std::string> NeighborsRecursive(const std::string_view pattern, int d);
std::set<std::string> ImmediateNeighbors(const std::string_view pattern);
std::set<std::string> NeighborsIterative(const std::string_view pattern, int d);

BIOUTILS_END_SUB_NAMESPACE(algorithms)

#endif // LIB_PATTERN_H