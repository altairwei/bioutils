#include <iostream>
#include <cstring>

using namespace std;

int PatternCount_1(char const *, char const*);
int PatternCount_2(char const *, char const*);

int main( int argc,      // Number of strings in array argv
          char *argv[],   // Array of command-line argument strings
          char *envp[] )  // Array of environment variable strings
{
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " TEXT PARTTERN" << endl;
        return 1;
    }

    cout << PatternCount_2(argv[1], argv[2]) << endl;
    return 0;
}

int PatternCount_1(char const *text, char const *parttern)
{
    int count = 0;
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

int PatternCount_2(char const *text, char const *parttern)
{
    int count = 0;
    long text_len = strlen(text);
    long parttern_len = strlen(parttern);

    for (int i = 0; i < text_len - parttern_len + 1; i++) {
        if (strncmp(&text[i], parttern, parttern_len) == 0) {
            count++;
        }
    }

    return count;
}