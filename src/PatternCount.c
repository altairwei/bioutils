#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>

#include "glib.h"

#define PROGRAM_NAME "wc"
#define ISNTP(C) ((C) == 'A' || (C) == 'T' || (C) == 'C' || (C) == 'G' \
        || (C) == 'a' || (C) == 't' || (C) == 'c' || (C) == 'g')

typedef unsigned long long hash_t;

unsigned int PatternCount_1(char *, char *);
unsigned int PatternCount_2(char *, char *);
unsigned int PatternCount_3(char *, char *);
hash_t hash_kmer(char *, size_t);
int bpton(char);
char *read_file(char *);
char *read_stdin();

static gint ALG = 2;

static GOptionEntry entries[] =
{
  { "algorithm", 'g', 0, G_OPTION_ARG_INT, &ALG, "Algorithm to be applied.", "N" },
  { NULL }
};

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
    GError *error = NULL;
    GOptionContext *context;

    context = g_option_context_new("PARTTERN [FILE]");
    g_option_context_add_main_entries(context, entries, NULL);
    if (!g_option_context_parse(context, &argc, &argv, &error)) {
        g_print("option parsing failed: %s\n", error->message);
        exit (1);
    }

    // Parse first positional argument.
    char *parttern = NULL;
    if (--argc < 1) {
        g_print(g_option_context_get_help(context, true, NULL));
        exit(1);
    } else {
        parttern = *++argv;
    }

    // Parse second positional argument.
    char *text;
    if (--argc < 1) {
        text = read_stdin();
    } else {
        text = read_file(*++argv);
    }

    unsigned int count;
    switch (ALG)
    {
    case 1:
        count = PatternCount_1(text, parttern);
        break;
    case 2:
        count = PatternCount_2(text, parttern);
        break;
    case 3:
        count = PatternCount_3(text, parttern);
        break;
    default:
        count = PatternCount_2(text, parttern);
        break;
    }

    printf("%i\n", count);

    free(text);
    return 0;
}

unsigned int
PatternCount_1(char *text, char *parttern)
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

unsigned int
PatternCount_2(char *text, char *parttern)
{
    unsigned int count = 0;
    long text_len = strlen(text);
    long parttern_len = strlen(parttern);

    for (int i = 0; i < text_len - parttern_len + 1; i++) {
        if (strncmp(&text[i], parttern, parttern_len) == 0) {
            count++;
        }
    }

    return count;
}

/**
 * @brief KM algorithm.
 * 
 * @param text 
 * @param parttern The max length of pattern is 32, which can be hashed in to `long long` type.
 * @return unsigned int 
 */
unsigned int
PatternCount_3(char *text, char *parttern)
{
    unsigned int count = 0;
    size_t t_len = strlen(text);
    size_t p_len = strlen(parttern);
    hash_t pattern_hash = hash_kmer(parttern, p_len);
    hash_t kmer_hash = hash_kmer(text, p_len); /* hash value of first kmer */

    // Triming non NTP characters
    for (; t_len > 1 && !ISNTP(text[t_len -1]); t_len--)
        continue;

    for (int i = 1; i < t_len - p_len + 1; i++) {
        // If hash values are matched then k-mer and pattern are matched.
        if (kmer_hash == pattern_hash)
            count++;
        // Compute hash of next k-mer
        //TODO: 用按位与来将高位清零，然后左移，再加上新的末尾hash
        kmer_hash = kmer_hash - (bpton(text[i-1]) << 2*(p_len - 1));
        kmer_hash = (kmer_hash << 2) + (bpton(text[i+ p_len - 1]));
    }

    return count;
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
        g_error("unknown base.");
        break;
    }

    return val;
}

hash_t
hash_kmer(char *kmer, size_t len)
{
    hash_t hash = 0;

    while (len-- > 0) {
        int val = bpton(*kmer++);
        hash += val << 2*len;
    }

    return hash;
}

//TODO: 将这些读取函数更改一下。函数返回实际读入的字符数，然后接受Buffer的指针，直接将其指向新的动态分配数组

/**
 * @brief Read file content to memory.
 * 
 * @param filename 
 * @return char* The pointer to text read in.
 */
char *
read_file(char *filename)
{
    char * buffer = 0;
    long length;
    FILE * fp = fopen (filename, "rb");

    if (fp)
    {
        fseek (fp, 0, SEEK_END);
        length = ftell(fp);
        fseek (fp, 0, SEEK_SET);
        buffer = (char *) malloc(length);
        if (buffer)
        {
            fread (buffer, 1, length, fp);
        }
        fclose (fp);
    }

    return buffer;
}

/**
 * @brief Read whole stdin text stream into memory.
 * 
 * @return char* The pointer to text read in.
 */
char *
read_stdin()
{
    size_t cap = 4096;
    size_t len = 0; 

    char *buffer = (char *)malloc(cap * sizeof (char));
    int c;

    while ((c = fgetc(stdin)) != EOF)
        {
            buffer[len] = c;

            if (++len == cap)
                // Make the output buffer twice its current size
                buffer = (char *)realloc(buffer, (cap *= 2) * sizeof (char));
        }

    // Trim off any unused bytes from the buffer
    buffer = (char *)realloc(buffer, (len + 1) * sizeof (char));

    buffer[len] = '\0';

    return buffer;
}