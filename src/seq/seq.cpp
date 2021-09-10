#include <string>
#include <map>
#include <algorithm>
#include <bitset>

#include <CLI/CLI.hpp>

#include "dataio.h"
#include "pattern.h"

using namespace std;

#define PROGRAM_NAME "bioseq - A tool to manipulate sequences."

static const map<char, char> base_comp = {
    {'A', 'T'}, {'T', 'A'},
    {'C', 'G'}, {'G', 'C'}
};


string
ReverseComplement(const string &oriSeq) noexcept(false)
{
    size_t len = oriSeq.length();
    string revSeq;
    revSeq.reserve(len);

    // Get complementary
    for (const char &nucleotide : oriSeq) {
        revSeq.push_back(base_comp.at(nucleotide));
    }

    // Get reversed
    std::reverse(revSeq.begin(), revSeq.end());

    return revSeq;
}


int
main(int argc, char *argv[], char *envp[])
{
    CLI::App app{PROGRAM_NAME};
    app.require_subcommand(1);

    string file_name = "-";
    bool do_reverse_complement = false;
    bool do_hash = false;
    CLI::App* conv_subapp = app.add_subcommand("conv", "Convert the sequence");
    conv_subapp->add_flag("-c,--reverse-complement", do_reverse_complement, "Get reverse complement sequence.");
    conv_subapp->add_flag("-H,--hash", do_hash, "Get hash number of the sequence.");
    conv_subapp->add_option("file", file_name, "file contains sequences.");

    CLI11_PARSE(app, argc, argv);

    string seq = bioutils::IO::read_input(file_name);

    seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
    seq.erase(std::remove(seq.begin(), seq.end(), '\r'), seq.end());

    string output = seq;

    if (do_reverse_complement) {
        try {
            output = ReverseComplement(seq);
        } catch(std::out_of_range &) {
            std::cout << "Unknown base." << std::endl;
        }
    }

    if (do_hash && !output.empty()) {
        auto hash = bioutils::algorithms::PatternToNumber(output);
        std::cout << std::bitset<8*sizeof(hash)>(hash) << std::endl;

        for (int i = 0; i < (8*sizeof(hash)/2 - output.length()); i++) {
            std::cout << "  ";
        }

        for (auto c : output) {
            std::cout << ' ' << c;
        }

        std::cout << std::endl;

    } else {
        std::cout << output << std::endl;
    }
}