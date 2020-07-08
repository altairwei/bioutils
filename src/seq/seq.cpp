#include <string>
#include <map>
#include <algorithm>

#include <argparse/argparse.hpp>

#include "dataio.h"

using namespace std;

#define PROGRAM_NAME "bpseq"

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
    argparse::ArgumentParser program(PROGRAM_NAME);
    program.add_argument("-c", "--reverse-complement").help("Get reverse complement sequence.")
        .default_value(false)
        .implicit_value(true);
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

    string seq = bioutils::IO::read_input(fileName);
    if (program["--reverse-complement"] == true) {
        seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
        seq.erase(std::remove(seq.begin(), seq.end(), '\r'), seq.end());
        try {
            std::cout << ReverseComplement(seq) << std::endl;
        } catch(std::out_of_range &) {
            std::cout << "Unknown base." << std::endl;
        }
        
    }
}