#include <string>
#include <set>

#include "gtest/gtest.h"

#include "pattern.h"

using namespace bioutils::algorithms;

TEST(TestFrequentWords, HandleNormalPattern) {
    // http://bioinformaticsalgorithms.com/data/debugdatasets/replication/FrequentWordsProblem.pdf
    std::set<std::string> expected = {"ATG"};
    std::set<std::string> output;
    std::string seq = "ATGATGATG";
    FrequentWords(seq, 3, output);
    EXPECT_EQ(expected, output);

    expected.clear();
    output.clear();
    seq.clear();
    expected = {"CATG", "GCAT"};
    seq = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
    FrequentWords(seq, 4, output);
    EXPECT_EQ(expected, output);

    expected.clear();
    output.clear();
    seq = "CGGAAGCGAGATTCGCGTGGCGTGATTCCGGCGGGCGTGGAG"
            "AAGCGAGATTCATTCAAGCCGGGAGGCGTGGCGTGGCGTGGC"
            "GTGCGGATTCAAGCCGGCGGGCGTGATTCGAGCGGCGGATTC"
            "GAGATTCCGGGCGTGCGGGCGTGAAGCGCGTGGAGGAGGCGT"
            "GGCGTGCGGGAGGAGAAGCGAGAAGCCGGATTCAAGCAAGCA"
            "TTCCGGCGGGAGATTCGCGTGGAGGCGTGGAGGCGTGGAGGC"
            "GTGCGGCGGGAGATTCAAGCCGGATTCGCGTGGAGAAGCGAG"
            "AAGCGCGTGCGGAAGCGAGGAGGAGAAGCATTCGCGTGATTC"
            "CGGGAGATTCAAGCATTCGCGTGCGGCGGGAGATTCAAGCGA"
            "GGAGGCGTGAAGCAAGCAAGCAAGCGCGTGGCGTGCGGCGGG"
            "AGAAGCAAGCGCGTGATTCGAGCGGGCGTGCGGAAGCGAGCGG";
    expected = {"CGGCGGGAGATT", "CGGGAGATTCAA", 
                "CGTGCGGCGGGA", "CGTGGAGGCGTG",
                "CGTGGCGTGCGG", "GCGTGCGGCGGG",
                "GCGTGGAGGCGT", "GCGTGGCGTGCG", 
                "GGAGAAGCGAGA", "GGAGATTCAAGC",
                "GGCGGGAGATTC", "GGGAGATTCAAG",
                "GTGCGGCGGGAG", "TGCGGCGGGAGA"};
    FrequentWords(seq, 12, output);
    EXPECT_EQ(expected, output);

    expected.clear();
    output.clear();
    seq.clear();
    seq = "TGGTAGCGACGTTGGTCCCGCCGCTTGAGAATCTGGATGAACA"
                        "TAAGCTCCCACTTGGCTTATTCAGAGAACTGGTCAACACTT"
                        "GTCTCTCCCAGCCAGGTCTGACCACCGGGCAACTTTTAGAG"
                        "CACTATCGTGGTACAAATAATGCTGCCAC";
    expected = {"TGG"};
    FrequentWords(seq, 3, output);
    EXPECT_EQ(expected, output);

    expected.clear();
    output.clear();
    seq.clear();
    seq = "CAGTGGCAGATGACATTTTGCTGGTCGACTGGTTACAACAACG"
            "CCTGGGGCTTTTGAGCAACGAGACTTTTCAATGTTGCACCG"
            "TTTGCTGCATGATATTGAAAACAATATCACCAAATAAATAA"
            "CGCCTTAGTAAGTAGCTTTT";
    expected = {"TTTT"};
    FrequentWords(seq, 4, output);
    EXPECT_EQ(expected, output);

    expected.clear();
    output.clear();
    seq.clear();
    seq = "ATACAATTACAGTCTGGAACCGGATGAACTGGCCGCAGGTTAA"
            "CAACAGAGTTGCCAGGCACTGCCGCTGACCAGCAACAACAA"
            "CAATGACTTTGACGCGAAGGGGATGGCATGAGCGAACTGAT"
            "CGTCAGCCGTCAGCAACGAGTATTGTTGCTGACCCTTAACA"
            "ATCCCGCCGCACGTAATGCGCTAACTAATGCCCTGCTG";
    expected = {"AACAA"};
    FrequentWords(seq, 5, output);
    EXPECT_EQ(expected, output);

    expected.clear();
    output.clear();
    seq.clear();
    seq = "CCAGCGGGGGTTGATGCTCTGGGGGTCACAAGATTGCATTTTT"
            "ATGGGGTTGCAAAAATGTTTTTTACGGCAGATTCATTTAAA"
            "ATGCCCACTGGCTGGAGACATAGCCCGGATGCGCGTCTTTT"
            "ACAACGTATTGCGGGGTAAAATCGTAGATGTTTTAAAATAG"
            "GCGTAAC";
    expected = {"AAAAT", "GGGGT", "TTTTA"};
    FrequentWords(seq, 5, output);
    EXPECT_EQ(expected, output);
}

TEST(TestFrequentWords, HandleEmptyPattern) {
    std::set<std::string> expected = {};
    std::set<std::string> output;
    FrequentWords("", 3, output);
    EXPECT_EQ(expected, output);

    expected.clear();
    output.clear();
    expected = {};
    FrequentWords("ATCGTAGTCGCTAG", 0, output);
    EXPECT_EQ(expected, output);
}

TEST(TestFrequentWords, HandleOutOfBounds) {
    std::set<std::string> expected = {};
    std::set<std::string> output;

    // k-mer can not be negative.
    expected.clear();
    output.clear();
    expected = {};
    FrequentWords("ATCGTAGTCGCTAG", -4, output);
    EXPECT_EQ(expected, output);

    // k-mer equals to sequence length
    expected.clear();
    output.clear();
    expected = {"ATG"};
    FrequentWords("ATG", 3, output);
    EXPECT_EQ(expected, output);

    // k-mer is larger than sequence length
    expected.clear();
    output.clear();
    expected = {};
    FrequentWords("ATG", 4, output);
    EXPECT_EQ(expected, output);
}