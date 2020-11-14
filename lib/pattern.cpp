#include "pattern.h"

using namespace std;

namespace bioutils {

namespace algorithms {

/**
 * @brief Brute force algorithm by hand.
 * 
 * @param text 
 * @param pattern 
 * @return size_t
 */
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

/**
 * @brief Brute force algorithm.
 * 
 * @param text 
 * @param pattern 
 * @return size_t
 */
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
 * @brief KM algorithm.
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
    hash_t pattern_hash = hash_kmer(pattern, p_len);
    hash_t kmer_hash = hash_kmer(text, p_len); /* hash value of first kmer */

    hash_t mask = 0;
    mask = ~((~mask) << (2*(p_len - 1)));

    // Triming non NTP characters
    for (; t_len > 1 && !is_ntp(text[t_len -1]); t_len--)
        continue;

    for (int i = 1; i < t_len - p_len + 1; i++) {
        // If hash values are matched then k-mer and pattern are matched.
        if (kmer_hash == pattern_hash)
            count++;
        // Compute hash of next k-mer
        kmer_hash = ((kmer_hash & mask) << 2) + bpton(text[i+ p_len - 1]);
    }

    return count;
}

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

void
PatternIndex(const char *text, const char *pattern, std::vector<size_t> &output)
{
    if (strlen(text) == 0 || strlen(pattern) == 0)
        return;
    find_do(text, pattern,
        [&](const size_t i, const char *, const char *){
            output.push_back(i);
        }
    );
}

void find_do(const char *text, const char *pattern,
    std::function<void(const size_t, const char *, const char *)> callback)
{
    size_t text_len = strlen(text);
    size_t pattern_len = strlen(pattern);

    for (size_t i = 0; i < text_len - pattern_len + 1; i++) {
        if (strncmp(&text[i], pattern, pattern_len) == 0) {
            callback(i, text, pattern);
        }
    }
}

std::set<std::string>
FrequentWords(const std::string text, const int k)
{
    size_t t_len = text.length();
    if (t_len == 0)
        return std::set<std::string>();
    if (k <= 0)
        return std::set<std::string>();
    if (k > t_len)
        return std::set<std::string>();

    // Total count of k-mer
    size_t n_kmer = t_len - k + 1;

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
FrequentWordsFast(const std::string text, const int k)
{
    return std::set<std::string>();
}


std::map<std::string, size_t>
FrequencyTable(const std::string text, const int k)
{
    return std::map<std::string, size_t>();
}


std::string
MaxMap(std::map<std::string, size_t> &input_map)
{
    return std::string();
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

/**
 * @brief Convert DNA base to number.
 * 
 * @param base 
 * @return int 
 */
inline int
bpton(char base)
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
        throw std::logic_error("Unknown base.");
        break;
    }

    return val;
}

hash_t
hash_kmer(const char *kmer, size_t len)
{
    hash_t hash = 0;

    while (len-- > 0) {
        int val = bpton(*kmer++);
        hash += val << 2*len;
    }

    return hash;
}

} // algorithms

} // bioutils