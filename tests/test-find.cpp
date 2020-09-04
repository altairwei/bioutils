#include <string>
#include <set>

#include "gtest/gtest.h"

#include "pattern.h"

using namespace bioutils::algorithms;

TEST(TestFrequentWords, HandleNormalPattern) {
    std::set<std::string> expected = {"ATG"};
    std::set<std::string> output;
    FrequentWords("ATGATGATG", 3, output);
    EXPECT_EQ(expected, output);

    expected.clear();
    output.clear();
    expected = {"CATG", "GCAT"};
    FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, output);
    EXPECT_EQ(expected, output);

    expected.clear();
    output.clear();
    std::string seq = "CGGAAGCGAGATTCGCGTGGCGTGATTCCGGCGGGCGTGGAG"
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
    expected = {"CGGCGGGAGATT", "CGGGAGATTCAA", "CGTGCGGCGGGA", "CGTGGAGGCGTG", "CGTGGCGTGCGG",
                "GCGTGCGGCGGG", "GCGTGGAGGCGT", "GCGTGGCGTGCG", "GGAGAAGCGAGA", "GGAGATTCAAGC",
                "GGCGGGAGATTC", "GGGAGATTCAAG", "GTGCGGCGGGAG", "TGCGGCGGGAGA"};
    FrequentWords(seq, 12, output);
    EXPECT_EQ(expected, output);
}

TEST(TestFrequentWords, HandleEmptyPattern) {
    std::set<std::string> expected = {};
    std::set<std::string> output;
    FrequentWords("", 3, output);
    EXPECT_EQ(expected, output);
}