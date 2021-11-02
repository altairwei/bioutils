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
#include "global.h"

using namespace std;

#define PROGRAM_NAME "biofind - A tool to find pattern in sequence."

using namespace BIOUTILS_NAMESPACE;

size_t
count(const string &text, const string &pattern, const int algorithm = 2)
{
    size_t count;
    switch (algorithm)
    {
    case 1:
        count = algorithms::PatternCount(text.c_str(), pattern.c_str(),
            algorithms::AlgorithmEfficiency::Slow);
        break;
    case 2:
        count = algorithms::PatternCount(text.c_str(), pattern.c_str(),
            algorithms::AlgorithmEfficiency::Fast);
        break;
    case 3:
        count = algorithms::PatternCount(text.c_str(), pattern.c_str(),
            algorithms::AlgorithmEfficiency::Faster);
        break;
    default:
        count = algorithms::PatternCount(text.c_str(), pattern.c_str(),
            algorithms::AlgorithmEfficiency::Fast);
        break;
    }

    return count;
}


int
main( int argc, char *argv[], char *envp[] )
{
    CLI::App app{PROGRAM_NAME};

    string pattern;
    string file_name = "-";
    int kmer;
    int algorithm = 2;

    app.add_option("file", file_name, "file contains sequences.");
    app.require_subcommand(1);

    CLI::App* count_subapp = app.add_subcommand("count", "Count pattern in the sequence");
    count_subapp->fallthrough();
    count_subapp->add_option("-p,--pattern", pattern, "k-mer pattern to count.")->required();
    count_subapp->add_option("-g,--algorithm", algorithm, "Algorithm to be applied.");
    count_subapp->callback([&]() {
        // Parse second positional argument.
        string text;
        if (file_name == "-") {
            text = IO::read_stdin();
        } else {
            text = IO::read_file(file_name);
        }

        cout << count(text, pattern, algorithm) << endl;
    });

    CLI::App* index_subapp = app.add_subcommand("index", "Get Index of Pattern in the Sequence");
    index_subapp->fallthrough();
    index_subapp->add_option("-p,--pattern", pattern, "k-mer pattern to index.")->required();
    index_subapp->callback([&]() {
        string text;
        if (file_name == "-") {
            text = IO::read_stdin();
        } else {
            text = IO::read_file(file_name);
        }

        std::vector<size_t> output = algorithms::PatternIndex(text.c_str(), pattern.c_str());

        for (size_t i : output) {
            cout << i << " ";
        }

        cout << endl;
    });


    CLI::App* freq_subapp = app.add_subcommand("freq", "Find Most Frequent k-mer");
    freq_subapp->fallthrough();
    freq_subapp->add_option("-k,--kmer", kmer, "Length of k-mer to find.")->required();
    freq_subapp->callback([&]() {
        string text;
        if (file_name == "-") {
            text = IO::read_stdin();
        } else {
            text = IO::read_file(file_name);
        }

        set<string> results = algorithms::FrequentWords(text, kmer);
        for (auto kmer : results) {
            cout << kmer << endl;
        }
    });

    int k, window_length, times;
    CLI::App* clumps_subapp = app.add_subcommand("clumps", "Find Patterns Forming Clumps in a String");
    clumps_subapp->fallthrough();
    clumps_subapp->add_option("-k,--k-mer", k, "a string pattern of length k")->required();
    clumps_subapp->add_option("-L,--window-length", window_length, "The length of a short interval of the genome")->required();
    clumps_subapp->add_option("-t,--times", times, "Pattern appears at least times")->required();
    clumps_subapp->callback([&]() {
        string seq = bioutils::IO::read_input(file_name);
        seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
        seq.erase(std::remove(seq.begin(), seq.end(), '\r'), seq.end());

        auto clumps = algorithms::FindClumps(seq, k, window_length, times);
        for (auto clp : clumps)
            std::cout << clp << " ";
        std::cout << std::endl;
    });

    CLI::App* skew_subapp = app.add_subcommand("skew", "Find a Position in a Genome Minimizing the Skew");
    skew_subapp->fallthrough();
    skew_subapp->callback([&] {
        string seq = bioutils::IO::read_input(file_name);
        seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
        seq.erase(std::remove(seq.begin(), seq.end(), '\r'), seq.end());

        auto locations = algorithms::FindMinimumSkew(seq);
        for (auto loc : locations)
            std::cout << loc << ' ';
        std::cout << std::endl;
    });

    CLI11_PARSE(app, argc, argv);

    return 0;
}