#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <streambuf>


std::string
read_file(const std::string &fileName)
{
    std::ifstream file(fileName);
    std::string content((std::istreambuf_iterator<char>(file)),
                    std::istreambuf_iterator<char>());
    return content;
}


/**
 * @brief Read file content to memory.
 * 
 * @param filename 
 * @return char* The pointer to text read in.
 */
char *
read_file(const char *filename)
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

std::string read_input(const std::string &argument)
{
    if (argument == "-") {
        return read_stdin();
    } else {
        return read_file(argument);
    }
}