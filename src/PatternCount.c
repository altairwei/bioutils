#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

unsigned int PatternCount_1(char const *, char const*);
unsigned int PatternCount_2(char const *, char const*);
char *read_file(char *);
char *read_stdin();

int
main( int argc, char *argv[], char *envp[] )
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s [OPTIONS] PARTTERN [FILE]\n", argv[0]);
        exit(1);
    }

    char *p_args[2] = {0, 0};
    int pn = 0;
    int ALG;
    char *cur_opt;
    while (--argc > 0) {
        cur_opt = *++argv;
        if (*cur_opt == '-') {
            // Parse options
            switch (*++cur_opt)
            {
            case 'a':
                // Choose algorithm
                --argc;
                ALG = atoi(*++argv);
                break;
            default:
                fprintf(stderr, "Unkown options: -%s", cur_opt);
                break;
            }
        } else {
            // Parse arguments: PARTTERN FILE
            if (pn < 2) {
                p_args[pn] = cur_opt;
            } else {
                fprintf(stderr, "too many arguments.\n");
                fprintf(stderr, "Usage: %s [OPTIONS] PARTTERN [FILE]\n", argv[0]);
                exit(1);
            }

            pn++;
        }

    }

    /* Read file or stdin to memory. */
    char *text;
    char *parttern = p_args[0];
    char *filename = p_args[1];
    if (filename) {
        text = read_file(filename);
    } else {
        text = read_stdin();
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
    default:
        count = PatternCount_2(text, parttern);
        break;
    }

    printf("%i\n", count);

    free(text);
    return 0;
}

unsigned int
PatternCount_1(char const *text, char const *parttern)
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
PatternCount_2(char const *text, char const *parttern)
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