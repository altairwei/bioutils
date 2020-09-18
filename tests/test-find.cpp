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

    expected.clear();
    output.clear();
    expected = {1, 3, 9};
    PatternIndex("GATATATGCATATACTT", "ATAT", output);
    EXPECT_EQ(expected, output);
}


TEST(TestPatternIndex, LongSequence) {
    std::vector<size_t> output;
    std::vector<size_t> expected = {
        19, 24, 38, 49, 56, 80, 128, 164, 186, 225, 230, 239, 387,
        403, 413, 419, 426, 471, 482, 508, 520, 604, 613, 618, 623,
        646, 651, 679, 684, 691, 713, 727, 747, 770, 777, 784, 801,
        829, 836, 841, 897, 947, 986, 991, 1011, 1036, 1075, 1148,
        1153, 1158, 1173, 1186, 1194, 1199, 1220, 1232, 1262, 1267,
        1303, 1329, 1369, 1386, 1395, 1407, 1444, 1467, 1472, 1477,
        1516, 1521, 1530, 1555, 1560, 1599, 1604, 1625, 1640, 1648,
        1653, 1666, 1680, 1698, 1728, 1733, 1745, 1770, 1800, 1805,
        1812, 1817, 1822, 1856, 1872, 1877, 1889, 1933, 1942, 1947,
        1952, 1972, 1983, 2004, 2016, 2021, 2032, 2041, 2046, 2073,
        2131, 2153, 2172, 2218, 2223, 2229, 2234, 2272, 2290, 2312,
        2430, 2440, 2460, 2465, 2486, 2497, 2547, 2560, 2595, 2645,
        2678, 2716, 2721, 2745, 2751, 2772, 2788, 2793, 2831, 2849,
        2854, 2860, 2865, 2900, 2905, 2911, 2916, 2941, 2947, 2960,
        2975, 2980, 2991, 2996, 3001, 3040, 3063, 3081, 3102, 3107,
        3112, 3124, 3129, 3142, 3152, 3157, 3188, 3193, 3216, 3224,
        3279, 3284, 3305, 3310, 3315, 3320, 3345, 3357, 3362, 3385,
        3397, 3402, 3418, 3431, 3445, 3517, 3526, 3537, 3580, 3585,
        3643, 3675, 3694, 3712, 3728, 3739, 3753, 3772, 3777, 3792,
        3797, 3824, 3835, 3847, 3852, 3857, 3862, 3877, 3882, 3888,
        3893, 3900, 3919, 3930, 3935, 3950, 4032, 4053, 4088
    };
    PatternIndex(
        "CCGAACACCCGTACACCGAACACCACACCACACCTTGCACACCACACCTACACC"
        "ACACACCACACCGGACACCCACACCCACACCACGAACACCGAGAGTACACCTAC"
        "ACCTGACACCGGGGATCGTCACACCAAGTGGTGATACACCCACACCCTTTACAC"
        "CTACACCACACCCGTACACCCTGAACACCACACCTAGAGAGTTGCACACCTCAC"
        "ACCGAAGGCACACCACACCATCCACACCATAAACACCGTTAACACCGTAGAACA"
        "CCCAGCACACCCTTACCGCATACACCGACGTTAGACACCCACACCGGCAGTCAC"
        "ACCGTACACCCATTCGGTCCACACCCTACACCGCCTGCCACACCTACTGAGTTA"
        "CACCGCATGACACCATTATCCGAACACACCAATATACACCAACACCATACACCA"
        "TTTAACACCCCAAAACACCGACACCGACACCGCAAGCCCACACCACACCCACAC"
        "CACAGACACCTACACCGTTTAGACACCAACACCGACACCACACCCCACACCCAA"
        "GACACCGCTACACCCTGCTGGACACCGACACCTACACCTCACACCGGACACCGC"
        "ACACACCGCCACACCAATCACACCACACCACACCAGTACAACACCGACACCTAC"
        "ACCACACCACACCCAGATACACCCACACCGGACACCACACCAAACACCATTACA"
        "CCCACACCGGTACACCACACCTCGTACACCAAGTAGACACCCAACACACCACAC"
        "CTTGATGACACCTGACACCATACACCAAACACCACACCGAGGTAGACACCACAC"
        "CGCCATCGACCACACCCTGACACCATACACCACACCACACCTAGTCGACACCCA"
        "CACCCTCACACCTGACACCCGCGGCATACACCCACACCACTTACACCTACACCG"
        "GGGGAAACACCGAAACACCTCAACACCGGACACCACACCTAAGACACCGGGCGA"
        "TACACCTGACCCTGACACCACACCACACCCAACACCCGAACACCACACCCAAAC"
        "CTTGACACCCACACCAAAACACCCTTTATTAAAACACCCCGACCACCAAACACC"
        "ACACCCCACACCGAACACCCACACCGCATACACCGGTCACACCTTATCTCGCCC"
        "ACACCCTACACCCCACACCACACCACACCACACCGTACCACACCACACCCCCAC"
        "ACCAAAACACCACACCACACCGGTTACACCCCACACCAACACCCACACCATTAC"
        "ACCTACACCGCAACACCTGCACACCACACCAAGACTGGAGACACCTACCACACC"
        "CTCGTTTACACCACCTGACACCTTACACCTCCGACACCAAAAACCCGTTGGGTC"
        "ATCGGATCAGGACACCTTTACACCACACCTTCGAGGACACCACGGACACCACAC"
        "CCCACACCACACCGGTACACCGCGTTCACACCTCACACCGACACCACACCCCCT"
        "GAACTGTATACACCACACCACACCAACCCAACACCCTAGAAGACACCTGCCACA"
        "CCTTACACCACACCACCGACACCAACACCCAAACACCTTTGACACACCACACCA"
        "ACACCGTACACCGCAACACCCGCATTACACCTTACACCACACCACACCCCCCTA"
        "CACCCACACCACACCCTCGGACACCAGTACACCACACCACAGATAGACACCATA"
        "CACCTTACACCACATACACCTTTCACACCACACCCACACCCCGCTTAGACACCG"
        "ACACCACACCACACCTGACACCACACCTCGCACACCGCCCTTACACCACACCCC"
        "AGCAGAAAACGAACACCCACACCACACCACACACCACACCACACCACACCGACA"
        "CCTGACACCTAAACACCCCCACACCACACCTCTCCAACACCACACCAACACCTA"
        "CACCAGAAAGACACCGACACCCGACACCCGCTGTTGTACACCCACACCATCGAC"
        "ACCACACCACACCACACCCTACACCGGCACACCATGCAAACACCACACACCTGG"
        "ACACCCACACCACACCGCACACCACACCACACCTACACCACCGACACCACACCA"
        "CACACCTACTCCACAACACCTACACCAAACACCCTACACCTACACCTACACCTA"
        "CATACACCTACACCTAATATTATGGACACCACACCTTCAGACACCGTACACCAC"
        "ACACCCTATGTTACACCACAGGCAGAATTTGACACCTCACACCCACACCCACAC"
        "CCGCACACCACACCAACACCACACCACACCCCCAACACCGCTCTTACACCTTAC"
        "ACCGACACCAACACCGACACCGACACCACACCCCAATATCCCTCACACCACACC"
        "TAACCAGTATACACCGTTGACAACACCCCAATTTACACCCCATACACCTCAGAC"
        "CACACACCGGACGGGCAACACCTACACCGATGTTACTTTACACCGGGCTCGCGG"
        "ACACCACTCGACACCAACACCCGACACCTTACACCACACCAGCTGCGTGAACAC"
        "CTACACCATCCCAACACCACACCGACACCGTATGGACACCTACACCTCGAGAGT"
        "TCCGCTAGAACACCACACCCATACACCATACACCGCGTACACCGAACACCGACA"
        "CCCACACCACACCCAATGACACCGATGACACCGGCTCGATACACCTACACCGAA"
        "CACCATCAGACACCGCGTACACCCAACACCTGACACCAACACCGCGGCACACCT"
        "AGTGACACCTACACCTACACCACACCATACACCCTACACCGATGAACACCAACA"
        "CCACTCTAAACACCCAGGACACCAACACACCTAGACACCACACCAACGACAGAG"
        "ACACCCTACACCTGCCAAGCTTTACACCATTGGTGAATCACACACCACACCAAC"
        "ACCACACCACACCGCTTACACCCGACCCGAAAACACCCACACCACACCAACACC"
        "ACACCACATTACTCCCGTTACACCTACACCAACACCACACCTTTACACCACACC"
        "CAGCAACACCACACCAAATGGACACCACACCACACCACACCTTAGCCGATGTGC"
        "CGACACCGCTGTCGTCACACCAGTGACACCTTAGCGTACACACCACACCCAACA"
        "CCTACACCACACCCGAAACACCTGACACCACACCACACCACACCCTACACCACA"
        "CCATGACCACACACCAGCCGACACCACACCATACACCTACACCGAAACACCTTT"
        "CTACACCACACCACACCTGAACACCTAGTCACACCACGACACCAACACCTGACC"
        "ACACCGGGGGACACCTTTGGAACGACACCTAACACCGCCACACCACACCACACC"
        "CGACACCTATAACACCACACCACACCACACCAAAGGCACACCTTAACACCCACA"
        "CCAAGGGCTACACCACACCACACCTCCAAAACAAGGGACACCACACCCAACACC"
        "ACACCACACCGCGTGGACACCACACCTTGACACCAAATTGTGCACACCACACCT"
        "GCACACCTTAAGAACGACACCGTCAGTACACCGAAACCCTATGACACCTGGGAC"
        "ACCTGGCACACCAACTACACCACACCCACACCACACACCTGGACACCGTTTCGC"
        "GAGTGTGGGTTGCTTGACACCACACCACACCGCGGCCTTACACCGCACACCGTA"
        "AACACCGTTGACACCTCATTACTCGACACCACACCGCACACCCACACCCGACAC"
        "CGAACACCACACCTGGGCATACACACCACACCGTACACCTACACCACACCTGTG"
        "CTACACCAGGGGTACACCACACCTAGTACACCACACCGATACACCCACACCACA"
        "CCACACCCACCAACACCACACCATCAAGAACACCCTATACACCCACACCACACC"
        "TACACCACACCCTACACCACACCACACCACACCATCGACACCTACACCACACCA"
        "ACACCACACCAAACACCACACCCACACCCGGACACCACACCCACACCACACCAT"
        "AACACCTAACACCACACACCTACACCTACTCTGCTAAACACCCAACACCTCTAC"
        "ACCCTGCCGACACCGCGACACCGGCGACACCCTGTTACACCACACCTCACACCT"
        "TCGACACCAGCCAGAGACACCGGACACCGACACCCCGAACACCAACACACCCGA", "ACACCA", output);
    EXPECT_EQ(expected, output);
}

TEST(TestPatternIndex, DonotReverseComplements) {
    // This dataset checks if your code is written correctly but is also taking into account
    // reverse complements, which we are not yet doing. Even though the reverse complement of
    // “ACAC” (which is “GTGT”) occurs in Genome, we only want to count occurrences of “ACAC”
    // specifically, which only occurs at index 4.

    std::vector<size_t> output;
    std::vector<size_t> expected = {4};
    PatternIndex("TTTTACACTTTTTTGTGTAAAAA", "ACAC", output);
    EXPECT_EQ(expected, output);
}

TEST(TestPatternIndex, OffByOneErrorsOfStart) {
    // This dataset checks for off-by-one errors at the beginning of Genome. Notice that “AAA”
    // occurs at the very beginning of Genome, so if you were to miss the first kmer of Genome, your
    // code would output the following: 46 51 74

    std::vector<size_t> output;
    std::vector<size_t> expected = {0, 46, 51, 74};
    PatternIndex(
        "AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAAT"
        "TAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAG"
        "CGCACAGACAGATAATAATTACAGAGTACACAACATCCAT", "AAA", output);
    EXPECT_EQ(expected, output);
}

TEST(TestPatternIndex, OffByOneErrorsOfEnd) {
    // This dataset checks for off-by-one errors at the endof Genome. Notice that “TTT” occurs
    // at the very end of Genome, so if you were to miss the last kmer of Genome, your code would
    // output the following: 88 93 98

    std::vector<size_t> output;
    std::vector<size_t> expected = {88, 92, 98, 132};
    PatternIndex(
        "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTC"
        "ACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTT"
        "TATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", "TTT", output);
    EXPECT_EQ(expected, output);
}

TEST(TestPatternIndex, OverlapPattern) {
    // This test dataset checks if your code correctly handles cases where instances of Pattern
    // overlap in Genome. In this case, if you did not count overlaps, you would only find the first and
    // last instances of ATA (ATATATA and ATATATA). However, there is indeed a third
    // occurrence, where the other two overlap (ATATATA).

    std::vector<size_t> output;
    std::vector<size_t> expected = {0, 2, 4};
    PatternIndex("ATATATA", "ATA", output);
    EXPECT_EQ(expected, output);
}

TEST(TestPatternIndex, EmptyInput) {
    std::vector<size_t> output;
    std::vector<size_t> expected;
    PatternIndex("", "ATA", output);
    EXPECT_EQ(expected, output);

    output.clear();
    expected.clear();
    PatternIndex("ATATATA", "", output);
    EXPECT_EQ(expected, output);

    output.clear();
    expected.clear();
    PatternIndex("", "", output);
    EXPECT_EQ(expected, output);
}

TEST(TestPatternIndex, NonNucleotide) {
    std::vector<size_t> output;
    std::vector<size_t> expected = {0};
    PatternIndex("ABCDEFG", "ABC", output);
    EXPECT_EQ(expected, output);
}