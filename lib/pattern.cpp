#include "pattern.h"

#include <cstring>
#include <cmath>
#include <algorithm>
#include <limits>
#include <array>
#include <iostream>
#include <cassert>

#include "exceptions.h"
#include "utils.h"

using namespace std;
using namespace bioutils::utils;

BIOUTILS_BEGIN_SUB_NAMESPACE(algorithms)

#define SubstrCount(text_length, pattern_length)  (text_length) - (pattern_length) + 1

static inline bool isPatternValid(const size_t seq_len, const size_t pattern_len)
{
    return seq_len != 0 && pattern_len > 0 && pattern_len <= seq_len;
}

static const char INT_TO_BASE[4] = {'A', 'C', 'G', 'T'};

static const int BASE_TO_INT[256] = {
    REPEAT_LIST_N(-1, 60), REPEAT_LIST_N(-1, 5),
//  A,  B, C,  D,  E,  F, G,  H,  I,  J,  K,  L,  M,  N,  O,  P,  Q,  R,  S, T,  U,  V,  W,  X,  Y,  Z
    0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1,
    REPEAT_LIST_N(-1, 6),
//  a,  b, c,  d,  e,  f, g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s, t,  u,  v,  w,  x,  y,  z
    0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1,
    REPEAT_LIST_N(-1, 100), REPEAT_LIST_N(-1, 30), REPEAT_LIST_N(-1, 3)
};

/*!
    \brief Convert DNA base to number
    
    \param base DNA base
    \return int Integer of DNA base
 */
static inline int NucleobaseToInt(const char base)
{
    int val;
    switch (toupper(base))
    {
        case 'A':
            val = 0;
            break;
        case 'C':
            val = 1;
            break;
        case 'G':
            val = 2;
            break;
        case 'T':
            val = 3;
            break;
        default:
            throw UnknownNucleotideError(base);
    }

    return val;
}

static inline char IntToNucleobase(const int i)
{
    char base;
    switch (i)
    {
        case 0:
            base = 'A';
            break;
        case 1:
            base = 'C';
            break;
        case 2:
            base = 'G';
            break;
        case 3:
            base = 'T';
            break;
        default:
            throw std::runtime_error("Only integer 0~3 can be converted to nucleobase symbol");
            break;
    }

    return base;
}

/*!
    Find all occurrences in \a text where \a pattern appears as a
    substring and then execute \a callback on each match.
 */
void find_do(const std::string_view text, const std::string_view pattern,
    std::function<void(const size_t, const std::string_view, const std::string_view)> callback)
{
    size_t text_len = text.length();
    size_t pattern_len = pattern.length();

    for (size_t i = 0; i < SubstrCount(text_len, pattern_len); i++) {
        if (text.substr(i, pattern_len) == pattern) {
            callback(i, text, pattern);
        }
    }
}

/*!
    Compute the Number of Times a Pattern Appears in a Text (Brute force)

    A k-mer is a string of length k. This function computes the number
    of times that a k-mer \a pattern appears as a substring of \a text.

    This version of PatternCount is kind of brute force algorithm.
 */
size_t PatternCount_BF(const std::string_view text, const std::string_view pattern)
{
    unsigned int count = 0;

    find_do(text, pattern,
        [&](const size_t i, const std::string_view, const std::string_view) {
            count++;
        }
    );

    return count;
}

/*!
    Another brute force version of PatternCount.
 */
#if 0
static size_t PatternCount_BFH(const std::string_view text, const std::string_view pattern)
{
    unsigned int count = 0;
    char const *pText;
    char const *pKmer;
    char const *pPattern;
    bool end = false;

    for (pText = text; *pText; pText++) {
        // compare k-mer and pattern
        for (pKmer = pText, pPattern = pattern; *pPattern; pKmer++, pPattern++) {
            // do not loop the remaining k-mer
            if (!*pKmer) {
                end = true;
            }
            if (*pKmer != *pPattern) {
                break;
            }

        }
        if (!*pPattern) {
            // pPattern point to '\0' means kmer == pattern
            count++;
        }
        // do not loop the remaining k-mer
        if (end)
            break;

    }

    return count;
}
#endif

/*!
    Compute the Number of Times a Pattern Appears in a Text (Rabin-Karp)

    Compute the number of times a \a pattern appears in a \a text. This one is
    implemented with Rabin-Karp algorithm. The max length of \a pattern is 32,
    which can be hashed in to `long long` type.
 */
size_t PatternCount_RK(const std::string_view text, const std::string_view pattern)
{
    size_t t_len = text.length();
    size_t p_len = pattern.length();

    if (!isPatternValid(t_len, p_len))
        return 0;

    auto first_kmer = text.substr(0, p_len);

    // Check first k-mer
    hash_t pattern_hash = PatternToNumber(pattern);
    hash_t kmer_hash = PatternToNumber(first_kmer);
    unsigned int count = kmer_hash == pattern_hash;

    hash_t mask = 0;
    mask = ~((~mask) << (2*(p_len - 1)));

    int base;
    for (int i = 1; i < SubstrCount(t_len, p_len); i++) {
        base = BASE_TO_INT[text[i+ p_len - 1]];
        if (base == -1) throw UnknownNucleotideError(base);
        // Compute hash of next k-mer
        kmer_hash = ((kmer_hash & mask) << 2) | base;
        // If hash values are matched then k-mer and pattern are matched.
        if (kmer_hash == pattern_hash)
            count++;
    }

    return count;
}

/*!
    Find the number of times that a k-mer appears as a substring of text.
 */
size_t  PatternCount(const std::string_view text, const std::string_view pattern, AlgorithmEfficiency algo) noexcept(false)
{
    switch (algo)
    {
    case AlgorithmEfficiency::Slow:
    case AlgorithmEfficiency::Fast:
        return PatternCount_BF(text, pattern);
        break;
    case AlgorithmEfficiency::Faster:
    case AlgorithmEfficiency::Fastest:
        return PatternCount_RK(text, pattern);
        break;
    default:
        throw std::runtime_error("Unknown algorithms.");
        break;
    }
}

/*!
    Find all starting positions in \a text where \a pattern appears as a
    substring.
 */
std::vector<size_t> PatternIndex(const std::string_view text, const std::string_view pattern)
{
    if (!isPatternValid(text.length(), pattern.length()))
        return std::vector<size_t>();
    std::vector<size_t> output;
    find_do(text, pattern,
        [&](const size_t i, const std::string_view, const std::string_view){
            output.push_back(i);
        }
    );
    return output;
}

/*!
    \brief Find All Approximate Occurrences of a Pattern in a String

    We say that a k-mer \a pattern appears as a substring of \a text with at
    most \a d mismatches if there is some k-mer substring pattern2 of \a text
    having \a d or fewer mismatches with \a pattern, i.e., HammingDistance(pattern, pattern2) ≤ d.
 */
std::vector<size_t> PatternIndexApproximate(const std::string_view text, const std::string_view pattern, const size_t d)
{
    size_t text_len = text.length();
    size_t pattern_len = pattern.length();

    std::vector<size_t> output;
    for (size_t i = 0; i < SubstrCount(text_len, pattern_len); i++) {
        if (HammingDistance(text.substr(i, pattern_len), pattern) <= d) {
            output.push_back(i);
        }
    }

    return output;
}

/**
 * @brief Find the most frequent k-mers in a string.
 * 
 * @param text Text to search
 * @param k Length of k-mer
 * @param algo Choose an algorithms
 * @return std::set<std::string> All most frequent k-mers in text.
 */
std::set<std::string> FrequentWords(
    const std::string_view text, const int k, AlgorithmEfficiency algo /*= Slow*/)
{
    switch (algo)
    {
    case AlgorithmEfficiency::Slow:
        return FrequentWordsSlow(text, k);
    case AlgorithmEfficiency::Fast:
        return FrequentWordsByPerfectHash(text, k);
    case AlgorithmEfficiency::Faster:
        return FrequentWordsBySorting(text, k);
    case AlgorithmEfficiency::Fastest:
        return FrequentWordsByStdHash(text, k);
    default:
        break;
    }

    return std::set<std::string>();
}

std::set<std::string>
FrequentWordsSlow(const std::string_view text, const int k)
{
    size_t t_len = text.length();

    if (!isPatternValid(t_len, k))
        return std::set<std::string>();

    // Total count of k-mer
    size_t n_kmer = SubstrCount(t_len, k);

    // Array to store counts for each k-mer.
    vector<int> kmer_count(n_kmer);

    // Array to store current k-mer.
    string cur_kmer;
    cur_kmer.reserve(k); 

    size_t max_count = 0;
    // Calculate counts of each k-mer
    for (int i = 0; i < n_kmer; i++) {
        cur_kmer = text.substr(i, k);
        int c = PatternCount_BF(text, cur_kmer);
        kmer_count[i] = c;
        if (c > max_count)
            max_count = c;
    }

    // Put most frequent patterns together
    std::set<std::string> output;
    for (int i = 0; i < n_kmer; i++) {
        if (kmer_count[i] == max_count) {
            std::string_view par = text.substr(i, k);
            // Set will remove duplicates automaticlly.
            output.insert(std::string(par));
        }
    }

    return output;
}

/*!
    \brief Find the Most Frequent Words in a String
    
    This version of FrequentWords is implemented with FrequencyTable()
 */
std::set<std::string> FrequentWordsByStdHash(const std::string_view text, const int k)
{
    if (!isPatternValid(text.length(), k))
        return std::set<std::string>();

    auto kmer_freq_table = FrequencyTable(text, k);
    size_t max = MaxMap(kmer_freq_table);

    std::set<std::string> max_freq;
    for (auto p : kmer_freq_table) {
        if (p.second == max)
            max_freq.insert(p.first);
    }

    return max_freq;
}

/*!
    \brief Find the Most Frequent Words in a String
    
    This version of FrequentWords is implemented with FrequencyArray()
 */
std::set<std::string> FrequentWordsByPerfectHash(const std::string_view text, const int k)
{
    if (!isPatternValid(text.length(), k))
        return std::set<std::string>();

    auto freq_array = FrequencyArray(text, k);
    size_t max = MaxArray(freq_array);

    std::set<std::string> max_freq;

    // Compare to FrequentWordsByStdHash, here too much time wasted looping
    // empty areas of FrequencyArray.
    for (int i = 0; i < freq_array.size(); i++) {
        if (freq_array[i] == max) {
            max_freq.insert(
                NumberToPatternBitwise(i, k)
            );
        }
    }

    return max_freq;
}

/*!
    \brief Find the Most Frequent Words in a String
    
    Given a string \a text , list all its \a k -mer in the order they appear in
    \a text, and convert each \a k -mer into an integer using \c PatternToNumber
    to produce an array index. Then we sort index. Since identical \a k -mer
    clump together in the sorted array, frequent \a k -mer are the longest runs
    of identical pattern hash in sorted index.

    This method is not as efficient as FrequentWordsByStdHash, because FrequencyTable
    is implemented with std::unordered_map and is very efficient. The length of
    FrequencyTable result will be shorter than "index/count" array of
    FrequentWordsBySorting and then sorting index array also increases the
    running time.
 */
std::set<std::string> FrequentWordsBySorting(const std::string_view text, const int k)
{
    if (!isPatternValid(text.length(), k))
        return std::set<std::string>();

    size_t n_kmer = SubstrCount(text.length(), k);

    // `index` is the pattern hash of each k-mer in text.
    std::vector<hash_t> index(n_kmer, 0);
    for (size_t i = 0; i < n_kmer; i++)
        index[i] = PatternToNumberBitwise(text.substr(i, k));

    // Sorted with the default operator<
    std::sort(index.begin(), index.end());

    size_t max = 1;
    // Frequent k-mers are the longest runs of identical pattern hash in
    // sorted index.
    std::vector<size_t> count(n_kmer, 1);
    for (size_t i = 1; i < n_kmer; i++) {
        if (index[i] == index[i-1])
            count[i] = count[i-1] + 1;
        if (count[i] > max)
            max = count[i];
    }

    std::set<std::string> max_freq;
    for (int i = 0; i < n_kmer; i++) {
        if (count[i] == max) {
            max_freq.insert(
                NumberToPatternBitwise(index[i], k)
            );
        }
    }

    return max_freq;
}

/*!
    \brief Find the Most Frequent Words with Mismatches in a String

    See https://rosalind.info/problems/ba1i/
 */
std::set<std::string> FrequentWordsWithMismatches(const std::string_view text, const int k, const int d, bool rev_comp)
{
    if (!isPatternValid(text.length(), k))
        return std::set<std::string>();

    auto kmer_freq_table = FrequencyTableWithMismatches(text, k, d, rev_comp);
    size_t max = MaxMap(kmer_freq_table);

    std::set<std::string> max_freq;
    for (auto p : kmer_freq_table) {
        if (p.second == max)
            max_freq.insert(p.first);
    }

    return max_freq;
}

/*!
    \brief Find the Most Frequent Words with Mismatches in a String

    See https://rosalind.info/problems/ba1i/
 */
std::set<std::string> FrequentWordsWithMismatchesBySorting(
    const std::string_view text, const int k, const int d, bool rev_comp)
{
    if (!isPatternValid(text.length(), k))
        return set<string>();

    size_t n_kmer = SubstrCount(text.length(), k);

    // 最大的不同在于，我们将生成的 Neighborhoods 当做所有在原始序列中出现过的 k-mer
    vector<string> neighborhoods;
    for (size_t i = 0; i < n_kmer; i++) {
        auto kmer = text.substr(i, k);
        auto neighbors = NeighborsRecursive(kmer, d);
        neighborhoods.reserve(neighbors.size());
        for (auto it = neighbors.begin(); it != neighbors.end(); ) {
            neighborhoods.push_back(
                std::move(neighbors.extract(it++).value()));
        }

        if (rev_comp) {
            auto rev_neighbors = NeighborsRecursive(ReverseComplement(kmer), d);
            neighborhoods.reserve(neighborhoods.size() + rev_neighbors.size());
            for (auto it = rev_neighbors.begin(); it != rev_neighbors.end(); ) {
                neighborhoods.push_back(
                    std::move(rev_neighbors.extract(it++).value()));
            }
        }

    }

    vector<hash_t> index(neighborhoods.size(), 0);
    for (size_t i = 0; i < neighborhoods.size(); i++)
        index[i] = PatternToNumberBitwise(neighborhoods[i]);

    std::sort(index.begin(), index.end());

    size_t max = 1;
    vector<size_t> count(index.size(), 1);
    for (size_t i = 1; i < index.size(); i++) {
        if (index[i] == index[i-1])
            count[i] = count[i-1] + 1;
        if (count[i] > max)
            max = count[i];
    }

    std::set<std::string> max_freq;
    for (int i = 0; i < count.size(); i++) {
        if (count[i] == max) {
            max_freq.insert(
                NumberToPatternBitwise(index[i], k)
            );
        }
    }

    return max_freq;
}

/*!
    Make a table corresponding to counting the number of occurrences of every k-mer.

    We would slide a length \a k window along \a text, and if the current k-mer
    substring of \a text does not occur in the table, then we would create a new 
    entry for it. Otherwise, we would add \c 1 to the entry corresponding to the
    current k-mer substring of \a text. We call this table the frequency table 
    for \a text and \a k.
 */
std::unordered_map<std::string, uint> FrequencyTable(const std::string_view text, const int k)
{
    size_t t_len = text.length();
    size_t n_kmer = SubstrCount(t_len, k);

    if (!isPatternValid(t_len, k))
        return std::unordered_map<std::string, uint>();

    std::unordered_map<std::string, uint> output;
    for (int i = 0; i < n_kmer; i++) {
        // Performing an insertion if such key does not already exist and
        // the mapped value is value-initialized (in this case, it's zero).
        output[std::string(text.substr(i, k))]++;
    }

    return output;
}

unordered_map<string, uint> FrequencyTableWithMismatches(
    const string_view text, const int k, const int d, bool rev_comp)
{
    size_t t_len = text.length();
    size_t n_kmer = SubstrCount(t_len, k);

    if (!isPatternValid(t_len, k))
        return std::unordered_map<std::string, uint>();

    std::unordered_map<std::string, uint> output;
    for (int i = 0; i < n_kmer; i++) {
        // Neighborhood relationships are mutual. When one k-mer shows up,
        // it means all the neighbors show up. Therefor, we can increase the
        // count for all of them.
        auto kmer = text.substr(i, k);
        auto neighbors = NeighborsRecursive(kmer, d);
        for (const auto &neigh : neighbors)
            output[neigh]++;

        if (rev_comp) {
            // When one k-mer shows up, it means the reverse complement of that
            // k-mer can be found on reverse complement of \a text. If k-mer
            // equals to its reverse complement, counts are doubled.
            auto rev_neighbors = NeighborsRecursive(ReverseComplement(kmer), d);
            for (const auto &rev : rev_neighbors)
                output[rev]++;
        }
    }

    return output;
}

/*!
    Generate the frequency array of a DNA string.

    Given an integer \a k, we define the frequency array of a string \a text as
    an array of length 4^k, where the i-th element of the array holds the
    number of times that the i-th k-mer (in the lexicographic order) appears in
    \a text . Therefore, the index of frequency array is hash value of a k-mer.
 */
std::vector<uint> FrequencyArray(const std::string_view text, const int k)
{
    size_t t_len = text.length();
    size_t n_kmer = SubstrCount(t_len, k);

    if (!isPatternValid(t_len, k))
        return std::vector<uint>();

    // Create a vector of size 4^k with all values as zero.
    std::vector<uint> freq_array(pow(4, k), 0);
 
    for (int i = 0; i < n_kmer; i++) {
        // Unlike std::map::operator[], this operator never inserts a new
        // element into the container. Accessing a nonexistent element
        // through this operator is undefined behavior.
        freq_array[PatternToNumber(text.substr(i, k))]++;
    }

    return freq_array;    
}

inline bool
is_ntp(char c)
{
  switch (c)
    {
    case 'A': case 'T': case 'C': case 'G':
    case 'a': case 't': case 'c': case 'g':
      return true;
    default:
      return false;
    }
}

const int MAX_HASHABLE_LENGTH = std::numeric_limits<hash_t>::digits / 2;
hash_t PatternToNumber(const std::string_view pattern, AlgorithmEfficiency algo)
{
    if (pattern.length() > MAX_HASHABLE_LENGTH)
        throw std::runtime_error(
            "The length of the pattern exceeds the maximum hashable length.");

    switch (algo)
    {
    case AlgorithmEfficiency::Slow:
        return PatternToNumberRecursive(pattern);
        break;
    case AlgorithmEfficiency::Fast:
    case AlgorithmEfficiency::Faster:
    case AlgorithmEfficiency::Fastest:
        return PatternToNumberBitwise(pattern);
        break;
    default:
        break;
    }

    return 0;
}

hash_t PatternToNumberBitwise(const std::string_view pattern)
{
    hash_t hash = 0;
    size_t len = pattern.length();
    auto p = pattern.begin();

    while (len-- > 0) {
        int val = BASE_TO_INT[*p++];
        hash |= val << 2*len;
    }

    return hash;
}

hash_t PatternToNumberRecursive(const std::string_view pattern)
{
    if (pattern.empty())
        return 0;

    char base = pattern.back();
    std::string_view prefix = pattern.substr(0, pattern.length() - 1);

    return 4 * PatternToNumberRecursive(prefix) + NucleobaseToInt(base);
}

std::string NumberToPatternBitwise(const hash_t number, const int length)
{
    std::string pattern;

    for (int len = length; len > 0; len--) {
        hash_t mask = 0;
        mask = ~((~mask) << (2*len));
        pattern.push_back(
            INT_TO_BASE[(number & mask) >> 2*(len - 1)]
        );
    }

    return pattern;
}

std::set<std::string> FindClumpsRaw(const std::string_view genome, int k, int window_length, int times)
{
    std::set<std::string> clumps;

    // Create a vector of size 4^k with all values as zero.
    std::vector<size_t> freq_array(pow(4, k), 0);
 
    // This is used to mark which pattern formed a clump. The pos of this
    // vector is hash value of a k-mer.
    std::vector<bool> is_clump(pow(4, k), false);
    for (size_t i = 0; i < SubstrCount(genome.length(), window_length); i++) {
        // Re-initialise to zero.
        std::fill(freq_array.begin(), freq_array.end(), 0);

        auto genome_window = genome.substr(i, window_length);
        for (int i = 0; i < SubstrCount(genome_window.length(), k); i++) {
            freq_array[PatternToNumber(genome_window.substr(i, k))]++;
        }

        for (size_t i = 0; i < freq_array.size(); i++) {
            if (freq_array[i] >= times)
                is_clump[i] = true;
        }
    }

    // Remove duplication
    for (size_t i = 0; i < is_clump.size(); i++) {
        if (is_clump[i])
            clumps.insert(NumberToPatternBitwise(i, k));
    }

    return clumps;
}

/*!
    We use FrequencyTable to implement FindClumps. It's faster than FrequencyArray.
 */
std::set<std::string> FindClumpsRaw2(const std::string_view genome, int k, int window_length, int times)
{
    std::set<std::string> clumps;

    // Loop each window in genome.
    for (size_t i = 0; i < SubstrCount(genome.length(), window_length); i++) {
        // Generate frequency table of each window.
        auto freq_table = FrequencyTable(genome.substr(i, window_length), k);
        // Find k-mers that appear a given times at least.
        for (auto p : freq_table) {
            if (p.second >= times) clumps.insert(p.first);
        }
    }

    return clumps;
}

/*!
    The max \a k is 32, which can be hashed in to \c hash_t type. The
    argument \a k is also limited by the available memory.
 */
std::set<std::string> FindClumpsBetterWithPerfectHash(const std::string_view genome, int k, int window_length, int times)
{
    std::set<std::string> clumps;
    // This is used to mark which pattern formed a clump. The pos of this
    // vector is hash value of a k-mer.
    std::vector<bool> is_clump(pow(4, k), false);

    // Handle First window
    auto freq_array = FrequencyArray(genome.substr(0, window_length), k);
    for (size_t i = 0; i < freq_array.size(); i++) {
        if (freq_array[i] >= times)
            is_clump[i] = true;
    }

    // Update following window
    for (size_t i = 1; i < SubstrCount(genome.length(), window_length); i++) {
        // Prior k-mer have been passed so we reduce it's frequency. Because
        // it's frequency is lower than last window, there is no need to check
        // it's number of occurrence.
        auto prior_kmer_hash = PatternToNumber(genome.substr(i-1, k));
        --freq_array[prior_kmer_hash];

        // Next k-mer is new so we increase it's frequncy, and check for it's
        // number of occurrence.
        auto next_kmer_hash = PatternToNumber(genome.substr(i + (window_length - k), k));
        if (++freq_array[next_kmer_hash] >= times && !is_clump[next_kmer_hash])
            is_clump[next_kmer_hash] = true;
    }

    for (size_t i = 0; i < is_clump.size(); i++) {
        if (is_clump[i])
            clumps.insert(NumberToPatternBitwise(i, k));
    }

    return clumps;
}

/*!
    By benchmark testing, this function is not as efficient on large data sets
    as FindClumpsBetterWithPerfectHash
 */
std::set<std::string> FindClumpsBetterWithStdHash(const std::string_view genome, int k, int window_length, int times)
{
    std::set<std::string> clumps;
    std::unordered_map<std::string, bool> is_clump;

    // Handle First window
    auto freq_table = FrequencyTable(genome.substr(0, window_length), k);
    for (auto p : freq_table)
        is_clump[p.first] = p.second >= times ? true : false;

    // Update following window
    for (size_t i = 1; i < SubstrCount(genome.length(), window_length); i++) {
        // FIXME in C++20: Heterogeneous lookup for unordered containers (transparent hashing)
        // see http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p0919r2.html
        auto prior_kmer = std::string(genome.substr(i-1, k));
        --freq_table[prior_kmer];

        auto next_kmer = std::string(genome.substr(i + (window_length - k), k));
        if (++freq_table[next_kmer] >= times && !is_clump[next_kmer])
            is_clump[next_kmer] = true;
    }

    for (auto p : is_clump)
        if (p.second) clumps.insert(p.first);

    return clumps;
}

/*!
    \brief Find patterns forming clumps in a string.

    Given integers \a window_length and \a times, a string pattern of length
    \a k forms an clump inside a (larger) string \a genome if there is an
    interval of \a genome of length \a window_length in which pattern appears
    at least \a times.

    We defined a \a k -mer as a "clump" if it appears many times within a short
    interval of the \a genome. We slide a window of fixed \a window_length
    along the \a genome, looking for a region where a \a k -mer appears given
    \a times in short succession.
 */
std::set<std::string> FindClumps(const std::string_view genome, int k, int window_length, int times, AlgorithmEfficiency algo)
{
    switch (algo)
    {
    case AlgorithmEfficiency::Slow:
        return FindClumpsRaw(genome, k, window_length, times);
    case AlgorithmEfficiency::Fast:
        return FindClumpsRaw2(genome, k, window_length, times);
    case AlgorithmEfficiency::Faster:
        return FindClumpsBetterWithPerfectHash(genome, k, window_length, times);
    case AlgorithmEfficiency::Fastest:
        return FindClumpsBetterWithStdHash(genome, k, window_length, times);
    default:
    {
        if (k > MAX_HASHABLE_LENGTH)
            return FindClumpsBetterWithStdHash(genome, k, window_length, times);
        else
            return FindClumpsBetterWithPerfectHash(genome, k, window_length, times);
    }
    }
}

/*!
    \brief Find a Position in a Genome Minimizing the Skew
    
    Define the skew of a DNA string \a genome as the difference between the
    total number of occurrences of 'G' and 'C' in \a genome. The skew should
    achieve a minimum at the position where the reverse half-strand ends and
    the forward half-strand begins.
 */
std::vector<size_t> FindMinimumSkew(const std::string_view genome)
{
    int accumulator = 0, min = 0;
    std::vector<int> skew(genome.length() + 1, 0);
    for (size_t i = 0; i < genome.length(); i++) {
        switch (genome[i])
        {
        case 'g':
        case 'G':
            skew[i+1] = ++accumulator;
            break;
        case 'c':
        case 'C':
            skew[i+1] = --accumulator;
            break;
        default:
            skew[i+1] = accumulator;
            break;
        }

        if (accumulator < min) {
            min = accumulator;
        }
    }

    std::vector<size_t> locations;
    for (size_t i = 0; i < skew.size(); i++) {
        if (skew[i] == min) locations.push_back(i);
    }

    return locations;
}

/*!
    \brief Compute the Hamming distance between two DNA strings

    We say that position i in k-mers p1 … pk and q1 … qk is a mismatch
    if pi ≠ qi. For example, CGAAT and CGGAC have two mismatches. The number of
    mismatches between strings p and q is called the Hamming distance between
    these strings and is denoted HammingDistance(p, q).
 */
size_t HammingDistance(const std::string_view pattern1, const std::string_view pattern2) noexcept(false)
{
    size_t d = 0;

    if (pattern1.length() != pattern2.length())
        throw std::runtime_error(
            "Can not compute Hamming distance between "
            "sequences of unequal length.");

    for (size_t i = 0; i < pattern1.length(); i++) {
        if (pattern1[i] != pattern2[i])
            d++;
    }

    return d;
}

const static char NUCLEOTIDES[4] = {'A', 'C', 'G', 'T'};
/*!
    \brief Generate the d-Neighborhood of a String

    The d-neighborhood Neighbors(Pattern, d) is the set of all k-mers whose
    Hamming distance from Pattern does not exceed \a d.
 */
std::set<std::string> NeighborsRecursive(const std::string_view pattern, int d)
{
    if (d == 0)
        return {std::string(pattern)};

    if (pattern.length() == 1)
        return {"A", "C", "G", "T"};

    std::set<std::string> neighborhood;
    auto suffix = pattern.substr(1, pattern.length() - 1);
    auto suffix_neighbors = NeighborsRecursive(suffix, d);
    for (auto text : suffix_neighbors) {
        // Hamming distance of text and suffix is either equal to d or less
        // than d. So we can genetate neighbors of raw text based on neighbors
        // of suffix.
        if (HammingDistance(suffix, text) < d) {
            for (auto n : NUCLEOTIDES)
                neighborhood.insert(n + text);
        } else {
            // Because the maximum allowed value of d has been reached, so we
            // just add first symbol of raw text to neighbors of suffix.
            neighborhood.insert(pattern.at(0) + text);
        }
    }

    return neighborhood;
}

/*!
    Generate the 1-neigborhood of \a pattern
 */
set<string> ImmediateNeighbors(const string_view pattern)
{
    set<string> neighborhood{string(pattern)};

    for (int i = 0; i < pattern.length(); i++) {
        const auto symbol = pattern.at(i);
        for (const auto N : NUCLEOTIDES) {
            if (symbol != N) {
                string neighbor(pattern);
                neighbor[i] = N;
                neighborhood.insert(neighbor);
            }
        }
    }

    return neighborhood;
}

/*!
    Iterative version of Neighbors avoiding recursive algorithm
 */
set<string> NeighborsIterative(const string_view pattern, int d)
{
    set<string> neighborhood{string(pattern)};

    for (int i = 1; i <= d; i++) {
        set<string> n_neighbors;
        for (const auto p : neighborhood) {
            auto nei = ImmediateNeighbors(p);
            n_neighbors.merge(nei);
        }
        neighborhood.merge(n_neighbors);
    }

    return neighborhood;
}

static const map<char, char> base_comp = {
    {'A', 'T'}, {'T', 'A'},
    {'C', 'G'}, {'G', 'C'}
};

string
ReverseComplement(const string_view oriSeq) noexcept(false)
{
    size_t len = oriSeq.length();
    string revSeq;
    revSeq.reserve(len);

    // Get complementary
    for (const char &nucleotide : oriSeq) {
        try {
            revSeq.push_back(base_comp.at(nucleotide));
        } catch(std::out_of_range &) {
            throw UnknownNucleotideError(nucleotide);
        }
    }

    // Get reversed
    std::reverse(revSeq.begin(), revSeq.end());

    return revSeq;
}

BIOUTILS_END_SUB_NAMESPACE(algorithms)