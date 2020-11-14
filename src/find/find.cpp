#include <cstring>
#include <cstdio>
#include <cctype>
#include <map>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>

#include <CLI/CLI.hpp>

#include "dataio.h"
#include "pattern.h"

using namespace std;

#define PROGRAM_NAME "biofind - A tool to find pattern in sequence."


size_t
count(const string &text, const string &pattern, const int algorithm = 2)
{
    size_t count;
    switch (algorithm)
    {
    case 1:
        count = bioutils::algorithms::PatternCount(text.c_str(), pattern.c_str(),
            bioutils::algorithms::PatternCountAlgorithms::BruteForceByHand);
        break;
    case 2:
        count = bioutils::algorithms::PatternCount(text.c_str(), pattern.c_str(),
            bioutils::algorithms::PatternCountAlgorithms::BruteForce);
        break;
    case 3:
        count = bioutils::algorithms::PatternCount(text.c_str(), pattern.c_str(),
            bioutils::algorithms::PatternCountAlgorithms::RabinKarp);
        break;
    default:
        count = bioutils::algorithms::PatternCount(text.c_str(), pattern.c_str(),
            bioutils::algorithms::PatternCountAlgorithms::BruteForce);
        break;
    }

    return count;
}


int
main( int argc, char *argv[], char *envp[] )
{
    CLI::App app{PROGRAM_NAME};
    app.require_subcommand(1);

    string pattern;
    string fileName = "-";
    int kmer;
    int algorithm = 2;

    CLI::App* count_subapp = app.add_subcommand("count", "Count pattern in the sequence");
    count_subapp->add_option("-p,--pattern", pattern, "k-mer pattern to count.")->required();
    count_subapp->add_option("-g,--algorithm", algorithm, "Algorithm to be applied.");
    count_subapp->add_option("file", fileName, "file contains sequences.");
    count_subapp->callback([&]() {
        // Parse second positional argument.
        string text;
        if (fileName == "-") {
            text = bioutils::IO::read_stdin();
        } else {
            text = bioutils::IO::read_file(fileName);
        }

        cout << count(text, pattern, algorithm) << endl;
    });

    CLI::App* index_subapp = app.add_subcommand("index", "Get index of pattern in the sequence");
    index_subapp->add_option("-p,--pattern", pattern, "k-mer pattern to index.")->required();
    index_subapp->add_option("file", fileName, "file contains sequences.");
    index_subapp->callback([&]() {
        string text;
        if (fileName == "-") {
            text = bioutils::IO::read_stdin();
        } else {
            text = bioutils::IO::read_file(fileName);
        }

        std::vector<size_t> output;
        bioutils::algorithms::PatternIndex(text.c_str(), pattern.c_str(), output);

        for (size_t i : output) {
            cout << i << " ";
        }

        cout << endl;
    });


    CLI::App* freq_subapp = app.add_subcommand("freq", "Find most frequent k-mer");
    freq_subapp->add_option("-k,--kmer", kmer, "Length of k-mer to find.")->required();
    freq_subapp->add_option("file", fileName, "file contains sequences.");
    freq_subapp->callback([&]() {
        string text;
        if (fileName == "-") {
            text = bioutils::IO::read_stdin();
        } else {
            text = bioutils::IO::read_file(fileName);
        }

        set<string> results = bioutils::algorithms::FrequentWords(text, kmer);
        for (auto kmer : results) {
            cout << kmer << endl;
        }
    });

    CLI11_PARSE(app, argc, argv);

    return 0;
}