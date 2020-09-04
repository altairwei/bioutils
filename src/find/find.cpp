#include <cstring>
#include <cstdio>
#include <cctype>
#include <map>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>

#include <argparse/argparse.hpp>

#include "dataio.h"
#include "pattern.h"

using namespace std;

#define PROGRAM_NAME "bpfind"

int
main( int argc, char *argv[], char *envp[] )
{
    argparse::ArgumentParser program(PROGRAM_NAME);
    program.add_argument("-g", "--algorithm").help("Algorithm to be applied.")
            .action([](const std::string &value) { return std::stoi(value); })
            .default_value(2);
    program.add_argument("-p", "--pattern").help("k-mer pattern to count.");
    program.add_argument("-k", "--kmer").help("Find most frequent k-mer.")
            .action([](const std::string &value) { return std::stoi(value); });
    program.add_argument("file").help("file contains sequences")
            .default_value(std::string("-"));

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cout << err.what() << std::endl;
        std::cout << program;
        exit(0);
    }

    string fileName = program.get<string>("file");
    //FIXME: use program.present instead when argparse v2.2 is aviable.
    try {
        string pattern = program.get<string>("--pattern");
        int ALG = program.get<int>("--algorithm");
        
        // Parse second positional argument.
        char *text;
        if (fileName == "-") {
            text = bioutils::IO::read_stdin();
        } else {
            text = bioutils::IO::read_file(fileName.c_str());
        }

        size_t count;
        switch (ALG)
        {
        case 1:
            count = bioutils::algorithms::PatternCount(text, pattern.c_str(),
                bioutils::algorithms::PatternCountAlgorithms::BruteForceByHand);
            break;
        case 2:
            count = bioutils::algorithms::PatternCount(text, pattern.c_str(),
                bioutils::algorithms::PatternCountAlgorithms::BruteForce);
            break;
        case 3:
            count = bioutils::algorithms::PatternCount(text, pattern.c_str(),
                bioutils::algorithms::PatternCountAlgorithms::RabinKarp);
            break;
        default:
            count = bioutils::algorithms::PatternCount(text, pattern.c_str(),
                bioutils::algorithms::PatternCountAlgorithms::BruteForce);
            break;
        }

        printf("%i\n", count);
        
        free(text);
    } catch (std::logic_error &) {
        // Option does not exist.
    }

    try {
        int k = program.get<int>("--kmer");

        string text;
        if (fileName == "-") {
            text = bioutils::IO::read_stdin();
        } else {
            text = bioutils::IO::read_file(fileName);
        }

        set<string> results;
        bioutils::algorithms::FrequentWords(text, k, results);

        for (auto kmer : results) {
            cout << kmer << endl;
        }

    } catch (std::logic_error &) {
        // Option does not exist.
    }


    return 0;
}