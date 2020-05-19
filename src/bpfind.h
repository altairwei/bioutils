#ifndef BPFIND_H
#define BPFIND_H

#include <map>
#include <string>
#include <vector>
#include <set>

#define ISNTP(C) ((C) == 'A' || (C) == 'T' || (C) == 'C' || (C) == 'G' \
        || (C) == 'a' || (C) == 't' || (C) == 'c' || (C) == 'g')

typedef unsigned long long hash_t;

unsigned int PatternCount_BFH(const char *, const char *);
unsigned int PatternCount_BF(const char *, const char *);
unsigned int PatternCount_KM(const char *, const char *);
void FrequentWords(const std::string text, const int k, std::set<std::string> &result);
hash_t hash_kmer(const char *, size_t);
bool is_ntp(char c);
int bpton(char);
char *read_file(char *);
char *read_stdin();

#endif //BPFIND_H