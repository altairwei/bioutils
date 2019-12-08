#include <iostream>
#include <string>

using namespace std;
int main( int argc,      // Number of strings in array argv
          char *argv[],   // Array of command-line argument strings
          char *envp[] )  // Array of environment variable strings
{
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " TEXT PARTTERN" << endl;
        return 1;
    }
}

int PatternCount(string text, string parttern)
{
    int count = 0;
    return count;
}