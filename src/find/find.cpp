#include <cstring>
#include <cstdio>
#include <cctype>
#include <map>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>

#include <argparse/argparse.hpp>

#include "find.h"
#include "dataio.h"

using namespace std;

#define PROGRAM_NAME "bpfind"

void emit_help()
{
    fprintf(stderr, "Usage: %s [OPTIONS] PARTTERN [FILE]\n", PROGRAM_NAME);
    exit(1);
}

void die(char *msg)
{
    fprintf(stderr, "Error: %s\n", msg);
    exit(1);
}

int
main( int argc, char *argv[], char *envp[] )
{
    argparse::ArgumentParser program(PROGRAM_NAME);
    program.add_argument("-g", "--algorithm").help("Algorithm to be applied.")
            .action([](const std::string &value) { return std::stoi(value); })
            .default_value(2);
    program.add_argument("-p", "--pattern").help("k-mer pattern to count.");
    program.add_argument("-k", "--kmer").help("Find most frequent k-mer.")
            .action([](const std::string &value) { return std::stoi(value); });
    program.add_argument("file").help("file contains sequences")
            .default_value(std::string("-"));

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cout << err.what() << std::endl;
        std::cout << program;
        exit(0);
    }

    string fileName = program.get<string>("file");
    //FIXME: use program.present instead when argparse v2.2 is aviable.
    try {
        string pattern = program.get<string>("--pattern");
        int ALG = program.get<int>("--algorithm");
        
        // Parse second positional argument.
        char *text;
        if (fileName == "-") {
            text = bioutils::IO::read_stdin();
        } else {
            text = bioutils::IO::read_file(fileName.c_str());
        }

        unsigned int count;
        switch (ALG)
        {
        case 1:
            count = PatternCount_BFH(text, pattern.c_str());
            break;
        case 2:
            count = PatternCount_BF(text, pattern.c_str());
            break;
        case 3:
            count = PatternCount_KM(text, pattern.c_str());
            break;
        default:
            count = PatternCount_BF(text, pattern.c_str());
            break;
        }

        printf("%i\n", count);
        
        free(text);
    } catch (std::logic_error &) {
        // Option does not exist.
    }

    try {
        int k = program.get<int>("--kmer");

        string text;
        if (fileName == "-") {
            text = bioutils::IO::read_stdin();
        } else {
            text = bioutils::IO::read_file(fileName);
        }

        set<string> results;
        FrequentWords(text, k, results);

        for (auto kmer : results) {
            cout << kmer << endl;
        }

    } catch (std::logic_error &) {
        // Option does not exist.
    }


    return 0;
}

/**
 * @brief Brute force algorithm by hand.
 * 
 * @param text 
 * @param parttern 
 * @return unsigned int 
 */
unsigned int
PatternCount_BFH(const char *text, const char *parttern)
{
    unsigned int count = 0;
    char const *pText;
    char const *pKmer;
    char const *pPattern;
    bool end = false;

    for (pText = text; *pText; pText++) {
        // compare k-mer and pattern
        for (pKmer = pText, pPattern = parttern; *pPattern; pKmer++, pPattern++) {
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
 * @param parttern 
 * @return unsigned int 
 */
unsigned int
PatternCount_BF(const char *text, const char *parttern)
{
    unsigned int count = 0;

    find_do(text, parttern,
        [&](const size_t i, const char *, const char *) {
            count++;
        }
    );

    return count;
}

void
PatternIndex(std::vector<size_t> indexes, const char *text, const char *parttern)
{
    find_do(text, parttern,
        [&](const size_t i, const char *, const char *){
            indexes.push_back(i);
        }
    );
}

void find_do(const char *text, const char *parttern,
    std::function<void(const size_t, const char *, const char *)> callback)
{
    size_t text_len = strlen(text);
    size_t parttern_len = strlen(parttern);

    for (size_t i = 0; i < text_len - parttern_len + 1; i++) {
        if (strncmp(&text[i], parttern, parttern_len) == 0) {
            callback(i, text, parttern);
        }
    }
}

/**
 * @brief KM algorithm.
 * 
 * @param text 
 * @param parttern The max length of pattern is 32, which can be hashed in to `long long` type.
 * @return unsigned int 
 */
unsigned int
PatternCount_KM(const char *text, const char *parttern)
{
    unsigned int count = 0;
    size_t t_len = strlen(text);
    size_t p_len = strlen(parttern);
    hash_t pattern_hash = hash_kmer(parttern, p_len);
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

void
FrequentWords(const std::string text, const int k, std::set<std::string> &result)
{
    size_t t_len = text.length();
    // Total count of k-mer
    // 究竟要如何理解 k-mer 的数量？
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
    for (int i = 0; i < n_kmer; i++) {
        if (kmer_count[i] == max_count) {
            string par = text.substr(i, k);
            // Set will remove duplicates automaticlly.
            result.insert(par);
        }
    }
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