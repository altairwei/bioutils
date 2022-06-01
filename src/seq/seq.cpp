#include <string>
#include <map>
#include <algorithm>
#include <bitset>

#include <CLI/CLI.hpp>

#include "dataio.h"
#include "pattern.h"
#include "exceptions.h"

using namespace std;
using namespace bioutils::utils;
using namespace bioutils;

#define PROGRAM_NAME "bioseq - A tool to manipulate sequences."

int
main(int argc, char *argv[], char *envp[])
{
    string file_name = "-";

    CLI::App app{PROGRAM_NAME};
    app.add_option("file", file_name, "file contains sequences.");
    app.require_subcommand(1);

    bool do_reverse_complement = false;
    bool do_hash = false;
    unsigned int hamming_distance = 0;
    CLI::App* conv_subapp = app.add_subcommand("conv", "Convert the sequence");
    conv_subapp->fallthrough();
    conv_subapp->add_flag("-c,--reverse-complement", do_reverse_complement, "Get reverse complement sequence.");
    conv_subapp->add_option("-d,--d-neighbors", hamming_distance, "Generate the d-Neighborhood of a String.");
    conv_subapp->add_flag("-H,--hash", do_hash, "Get hash number of the sequence.");
    conv_subapp->callback([&]() {
        string seq = bioutils::IO::read_input(file_name);

        seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
        seq.erase(std::remove(seq.begin(), seq.end(), '\r'), seq.end());

        string output = seq;

        if (do_reverse_complement)
            output = algorithms::ReverseComplement(seq);

        if (do_hash && !output.empty()) {
            auto hash = bioutils::algorithms::PatternToNumber(output);
            std::cout << hash << std::endl;
            std::cout << std::bitset<8*sizeof(hash)>(hash) << std::endl;

            for (int i = 0; i < (8*sizeof(hash)/2 - output.length()); i++) {
                std::cout << "  ";
            }

            for (auto c : output) {
                std::cout << ' ' << c;
            }

            std::cout << std::endl;

        } else if (hamming_distance > 0) {
            auto neighbors = bioutils::algorithms::NeighborsRecursive(seq, hamming_distance);
            for (auto n : neighbors)
                std::cout << n << std::endl;
        } else {
            std::cout << output << std::endl;
        }
    });

    int line_length = 50;
    std::string line_prefix;
    std::string line_suffix;
    CLI::App* format_subapp = app.add_subcommand("format", "Format the sequences with given style.");
    format_subapp->fallthrough();
    format_subapp->add_option("-L,--line-length", line_length, "How many characters per line.");
    format_subapp->add_option("-P,--line-prefix", line_prefix, "Add prefix to each line.");
    format_subapp->add_option("-S,--line-suffix", line_suffix, "Add suffix to each line.");
    format_subapp->callback([&]() {
        string seq = bioutils::IO::read_input(file_name);

        seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
        seq.erase(std::remove(seq.begin(), seq.end(), '\r'), seq.end());

        int count = 0;
        for (auto bp : seq) {
            if (count == 0 && !line_prefix.empty())
                std::cout << line_prefix;

            std::cout << bp;

            if (++count == line_length) {
                if (!line_suffix.empty())
                    std::cout << line_suffix;
                std::cout << std::endl;
                count = 0;
            }
        }

        if (count != 0) {
            if (!line_suffix.empty())
                std::cout << line_suffix;

            std::cout << std::endl;
        }
    });

    CLI11_PARSE(app, argc, argv);
}