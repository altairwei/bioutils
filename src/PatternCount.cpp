#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;

int PatternCount_1(char const *, char const*);
int PatternCount_2(char const *, char const*);

int
main( int argc, char *argv[], char *envp[] )
{
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " TEXT PARTTERN" << endl;
        return 1;
    }

    cout << PatternCount_2(argv[1], argv[2]) << endl;
    return 0;
}

int
PatternCount_1(char const *text, char const *parttern)
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

int
PatternCount_2(char const *text, char const *parttern)
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