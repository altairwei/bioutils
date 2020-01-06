#ifndef BPFIND_H
#define BPFIND_H

#define ISNTP(C) ((C) == 'A' || (C) == 'T' || (C) == 'C' || (C) == 'G' \
        || (C) == 'a' || (C) == 't' || (C) == 'c' || (C) == 'g')

typedef unsigned long long hash_t;

unsigned int PatternCount_1(const char *, const char *);
unsigned int PatternCount_2(const char *, const char *);
unsigned int PatternCount_3(const char *, const char *);
hash_t hash_kmer(const char *, size_t);
bool is_ntp(char c);
int bpton(char);
char *read_file(char *);
char *read_stdin();

#endif //BPFIND_H