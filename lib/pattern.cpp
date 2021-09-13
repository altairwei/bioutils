#include "pattern.h"

#include <cstring>

#include "exceptions.h"

using namespace std;
using namespace bioutils::utils;

BIOUTILS_BEGIN_SUB_NAMESPACE(algorithms)

#define PatternLoopCount(text_length, pattern_length)  (text_length) - (pattern_length) + 1

/** Brute force algorithm by hand. */
static
size_t
PatternCount_BFH(const char *text, const char *pattern)
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

/** Brute force algorithm. */
static
size_t
PatternCount_BF(const char *text, const char *pattern)
{
    unsigned int count = 0;

    find_do(text, pattern,
        [&](const size_t i, const char *, const char *) {
            count++;
        }
    );

    return count;
}

/**
 * @brief PatternCount implemented with RK algorithm.
 * 
 * @param text 
 * @param pattern The max length of pattern is 32, which can be hashed in to `long long` type.
 * @return size_t
 */
static
size_t
PatternCount_RK(const char *text, const char *pattern)
{
    unsigned int count = 0;
    size_t t_len = strlen(text);
    size_t p_len = strlen(pattern);
    hash_t pattern_hash = PatternToNumber(pattern);
    hash_t kmer_hash = PatternToNumber(text); /* hash value of first kmer */

    hash_t mask = 0;
    mask = ~((~mask) << (2*(p_len - 1)));

    // Triming non NTP characters
    for (; t_len > 1 && !is_ntp(text[t_len -1]); t_len--)
        continue;

    for (int i = 1; i < PatternLoopCount(t_len, p_len); i++) {
        // If hash values are matched then k-mer and pattern are matched.
        if (kmer_hash == pattern_hash)
            count++;
        // Compute hash of next k-mer
        kmer_hash = ((kmer_hash & mask) << 2) + NucleobaseToInt(text[i+ p_len - 1]);
    }

    return count;
}

/**
 * @brief Find the number of times that a k-mer appears as a substring of text.
 * 
 * @param text Text to find pattern.
 * @param pattern k-mer
 * @param algo 
 * @return size_t Number of times
 */
size_t 
PatternCount(const char *text, const char *pattern, enum PatternCountAlgorithms algo) noexcept(false)
{
    switch (algo)
    {
    case PatternCountAlgorithms::BruteForce:
        return PatternCount_BF(text, pattern);
        break;
    case PatternCountAlgorithms::BruteForceByHand:
        return PatternCount_BFH(text, pattern);
        break;
    case PatternCountAlgorithms::RabinKarp:
        return PatternCount_RK(text, pattern);
        break;
    default:
        throw std::runtime_error("Unknown algorithms.");
        break;
    }
}

static inline bool isPatternValid(const size_t seq_len, const size_t pattern_len)
{
    return seq_len != 0 && pattern_len > 0 && pattern_len <= seq_len;
}

std::vector<size_t>
PatternIndex(const char *text, const char *pattern)
{
    if (!isPatternValid(strlen(text), strlen(pattern)))
        return std::vector<size_t>();
    std::vector<size_t> output;
    find_do(text, pattern,
        [&](const size_t i, const char *, const char *){
            output.push_back(i);
        }
    );
    return output;
}

void find_do(const char *text, const char *pattern,
    std::function<void(const size_t, const char *, const char *)> callback)
{
    size_t text_len = strlen(text);
    size_t pattern_len = strlen(pattern);

    for (size_t i = 0; i < PatternLoopCount(text_len, pattern_len); i++) {
        if (strncmp(&text[i], pattern, pattern_len) == 0) {
            callback(i, text, pattern);
        }
    }
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
    const std::string &text, const int k, AlgorithmEfficiency algo /*= Slow*/)
{
    switch (algo)
    {
    case AlgorithmEfficiency::Slow:
        return FrequentWordsSlow(text, k);
        break;
    case AlgorithmEfficiency::Fast:
    case AlgorithmEfficiency::Faster:
    case AlgorithmEfficiency::Fastest:
        return FrequentWordsFast(text, k);
        break;
    default:
        break;
    }

    return std::set<std::string>();
}

std::set<std::string>
FrequentWordsSlow(const std::string &text, const int k)
{
    size_t t_len = text.length();

    if (!isPatternValid(t_len, k))
        return std::set<std::string>();

    // Total count of k-mer
    size_t n_kmer = PatternLoopCount(t_len, k);

    // Array to store counts for each k-mer.
    vector<int> kmer_count(n_kmer);

    // Array to store current k-mer.
    string cur_kmer;
    cur_kmer.reserve(k); 

    size_t max_count = 0;
    // Calculate counts of each k-mer
    for (int i = 0; i < n_kmer; i++) {
        cur_kmer = text.substr(i, k);
        int c = PatternCount_BF(text.c_str(), cur_kmer.c_str());
        kmer_count[i] = c;
        if (c > max_count)
            max_count = c;
    }

    // Put most frequent patterns together
    std::set<std::string> output;
    for (int i = 0; i < n_kmer; i++) {
        if (kmer_count[i] == max_count) {
            string par = text.substr(i, k);
            // Set will remove duplicates automaticlly.
            output.insert(par);
        }
    }

    return output;
}


std::set<std::string>
FrequentWordsFast(const std::string &text, const int k)
{
    size_t t_len = text.length();

    if (!isPatternValid(t_len, k))
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

/**
 * @brief Make a table corresponding to counting the number of occurrences of every k-mer.
 * 
 *  We would slide a length-k window Text, and if the current k-mer substring
 *  of text does not occur in the table, then we would create a new entry for
 *  it. Otherwise, we would add 1 to the entry corresponding to the current
 *  k-mer substring of Text. We call this table the frequency table for Text
 *  and k.
 * 
 * @param text Text to search
 * @param k Length of k-mer
 * @return std::map<std::string, size_t> A string-number map
 */
std::map<std::string, size_t>
FrequencyTable(const std::string &text, const int k)
{
    size_t t_len = text.length();
    size_t n_kmer = PatternLoopCount(t_len, k);

    if (!isPatternValid(t_len, k))
        return std::map<std::string, size_t>();

    std::map<std::string, size_t> output;
    for (int i = 0; i < n_kmer; i++) {
        // Performing an insertion if such key does not already exist and
        // the mapped value is value-initialized (in this case, it's zero).
        output[text.substr(i, k)]++;
    }

    return output;
}

/**
 * @brief Find the maximum value of a map
 * 
 * @param input_map Input map
 * @return size_t Maxium value of the map
 */
size_t
MaxMap(const std::map<std::string, size_t> &input_map) noexcept(false)
{
    if (input_map.empty())
        throw std::runtime_error("input map is empty.");

    auto p = input_map.cbegin();
    size_t max = (*p++).second;

    while (p != input_map.cend()) {
        size_t cur = (*p++).second;
        if (cur > max)
            max = cur;
    }

    return max;
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

/*!
    \brief Convert DNA base to number
    
    \param base DNA base
    \return int Integer of DNA base
 */
inline int NucleobaseToInt(char base)
{
    int val;
    switch (toupper(base))
    {
        case 'A':
            val = 0;
            break;
        case 'T':
            val = 1;
            break;
        case 'C':
            val = 2;
            break;
        case 'G':
            val = 3;
            break;
        default:
            throw UnknownNucleotideError(base);
    }

    return val;
}

hash_t
PatternToNumber(const std::string &kmer)
{
    hash_t hash = 0;
    size_t len = kmer.length();
    auto p = kmer.begin();

    while (len-- > 0) {
        int val = NucleobaseToInt(*p++);
        hash += val << 2*len;
    }

    return hash;
}



BIOUTILS_END_SUB_NAMESPACE(algorithms)