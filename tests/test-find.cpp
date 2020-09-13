#include <string>
#include <set>
#include <vector>

#include "gtest/gtest.h"

#include "pattern.h"

using namespace bioutils::algorithms;

TEST(TestFrequentWords, HandleNormalPattern) {
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
}

TEST(TestFrequentWords, CountFirstKmer) {
    // This dataset just checks if you’re counting the first kmer
    //  in Text(TGG in this example).
    std::string seq = "TGGTAGCGACGTTGGTCCCGCCGCTTGAGAATCTGGATGAACA"
        "TAAGCTCCCACTTGGCTTATTCAGAGAACTGGTCAACACTT"
        "GTCTCTCCCAGCCAGGTCTGACCACCGGGCAACTTTTAGAG"
        "CACTATCGTGGTACAAATAATGCTGCCAC";
    std::set<std::string> expected = {"TGG"};
    std::set<std::string> output;
    FrequentWords(seq, 3, output);
    EXPECT_EQ(expected, output);
}

TEST(TestFrequentWords, CountLastKmer) {
    // This dataset just checks if you’re counting the last kmer
    //  in Text (TTTT in this example).
    std::string seq = "CAGTGGCAGATGACATTTTGCTGGTCGACTGGTTACAACAACG"
        "CCTGGGGCTTTTGAGCAACGAGACTTTTCAATGTTGCACCG"
        "TTTGCTGCATGATATTGAAAACAATATCACCAAATAAATAA"
        "CGCCTTAGTAAGTAGCTTTT";
    std::set<std::string> expected = {"TTTT"};
    std::set<std::string> output;
    FrequentWords(seq, 4, output);
    EXPECT_EQ(expected, output);
}

TEST(TestFrequentWords, OverlappingOccurrences) {
    // This dataset checks if your code correctly handles cases
    //  where there are overlapping occurrences of Pattern
    //  throughout Text.
    std::string seq = "ATACAATTACAGTCTGGAACCGGATGAACTGGCCGCAGGTTAA"
        "CAACAGAGTTGCCAGGCACTGCCGCTGACCAGCAACAACAA"
        "CAATGACTTTGACGCGAAGGGGATGGCATGAGCGAACTGAT"
        "CGTCAGCCGTCAGCAACGAGTATTGTTGCTGACCCTTAACA"
        "ATCCCGCCGCACGTAATGCGCTAACTAATGCCCTGCTG";
    std::set<std::string> expected = {"AACAA"};
    std::set<std::string> output;
    FrequentWords(seq, 5, output);
    EXPECT_EQ(expected, output);
}

TEST(TestFrequentWords, OutputAllMostFrequentKmer) {
    // This test dataset checks if your code correctly handles ties
    //  (i.e. your code actually outputs ALL “most frequent” kmers,
    //  and not just a single “most frequent” kmer).
    std::string seq = "CCAGCGGGGGTTGATGCTCTGGGGGTCACAAGATTGCATTTTT"
        "ATGGGGTTGCAAAAATGTTTTTTACGGCAGATTCATTTAAA"
        "ATGCCCACTGGCTGGAGACATAGCCCGGATGCGCGTCTTTT"
        "ACAACGTATTGCGGGGTAAAATCGTAGATGTTTTAAAATAG"
        "GCGTAAC";
    std::set<std::string> expected = {"AAAAT", "GGGGT", "TTTTA"};
    std::set<std::string> output;
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


TEST(TestPatternIndex, HandleNormalInput) {
    std::vector<size_t> output;
    std::vector<size_t> expected;

    expected = {0, 3, 9};
    PatternIndex("ATGATGAAAATG", "ATG", output);
    EXPECT_EQ(expected, output);

    expected.clear();
    output.clear();
    expected = {0, 2, 12, 14};
    PatternIndex("ATATAAAAGGGGATATA", "ATA", output);
    EXPECT_EQ(expected, output);
}