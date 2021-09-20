#include <string>
#include <set>
#include <vector>

#include "gtest/gtest.h"

#include "pattern.h"

namespace {

using namespace bioutils::algorithms;
using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::ValuesIn;

typedef std::map<std::string, size_t> StrNumDict;

class TestFrequentWords : public TestWithParam<AlgorithmEfficiency> {
  // You can implement all the usual fixture class members here.
  // To access the test parameter, call GetParam() from class
  // TestWithParam<T>.
};

TEST_P(TestFrequentWords, HandleNormalPattern) {
    EXPECT_EQ(
        FrequentWords("ATGATGATG", 3, GetParam()),
        std::set<std::string>({"ATG"})
    );

    EXPECT_EQ(
        FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, GetParam()),
        std::set<std::string>({"CATG", "GCAT"})
    );


    EXPECT_EQ(
        FrequentWords(
            "CGGAAGCGAGATTCGCGTGGCGTGATTCCGGCGGGCGTGGAG"
            "AAGCGAGATTCATTCAAGCCGGGAGGCGTGGCGTGGCGTGGC"
            "GTGCGGATTCAAGCCGGCGGGCGTGATTCGAGCGGCGGATTC"
            "GAGATTCCGGGCGTGCGGGCGTGAAGCGCGTGGAGGAGGCGT"
            "GGCGTGCGGGAGGAGAAGCGAGAAGCCGGATTCAAGCAAGCA"
            "TTCCGGCGGGAGATTCGCGTGGAGGCGTGGAGGCGTGGAGGC"
            "GTGCGGCGGGAGATTCAAGCCGGATTCGCGTGGAGAAGCGAG"
            "AAGCGCGTGCGGAAGCGAGGAGGAGAAGCATTCGCGTGATTC"
            "CGGGAGATTCAAGCATTCGCGTGCGGCGGGAGATTCAAGCGA"
            "GGAGGCGTGAAGCAAGCAAGCAAGCGCGTGGCGTGCGGCGGG"
            "AGAAGCAAGCGCGTGATTCGAGCGGGCGTGCGGAAGCGAGCGG", 12, GetParam()),
        std::set<std::string>({
            "CGGCGGGAGATT", "CGGGAGATTCAA", 
            "CGTGCGGCGGGA", "CGTGGAGGCGTG",
            "CGTGGCGTGCGG", "GCGTGCGGCGGG",
            "GCGTGGAGGCGT", "GCGTGGCGTGCG", 
            "GGAGAAGCGAGA", "GGAGATTCAAGC",
            "GGCGGGAGATTC", "GGGAGATTCAAG",
            "GTGCGGCGGGAG", "TGCGGCGGGAGA"})
    );
}

TEST_P(TestFrequentWords, CountFirstKmer) {
    // This dataset just checks if you’re counting the first kmer
    //  in Text(TGG in this example).
    EXPECT_EQ(
        FrequentWords(
            "TGGTAGCGACGTTGGTCCCGCCGCTTGAGAATCTGGATGAACA"
            "TAAGCTCCCACTTGGCTTATTCAGAGAACTGGTCAACACTT"
            "GTCTCTCCCAGCCAGGTCTGACCACCGGGCAACTTTTAGAG"
            "CACTATCGTGGTACAAATAATGCTGCCAC", 3, GetParam()),
        std::set<std::string>({"TGG"})
    );
}

TEST_P(TestFrequentWords, CountLastKmer) {
    // This dataset just checks if you’re counting the last kmer
    //  in Text (TTTT in this example).
    EXPECT_EQ(
        FrequentWords(
            "CAGTGGCAGATGACATTTTGCTGGTCGACTGGTTACAACAACG"
            "CCTGGGGCTTTTGAGCAACGAGACTTTTCAATGTTGCACCG"
            "TTTGCTGCATGATATTGAAAACAATATCACCAAATAAATAA"
            "CGCCTTAGTAAGTAGCTTTT", 4, GetParam()),
        std::set<std::string>({"TTTT"})
    );
}

TEST_P(TestFrequentWords, OverlappingOccurrences) {
    // This dataset checks if your code correctly handles cases
    //  where there are overlapping occurrences of Pattern
    //  throughout Text.
    EXPECT_EQ(
        FrequentWords(
            "ATACAATTACAGTCTGGAACCGGATGAACTGGCCGCAGGTTAA"
            "CAACAGAGTTGCCAGGCACTGCCGCTGACCAGCAACAACAA"
            "CAATGACTTTGACGCGAAGGGGATGGCATGAGCGAACTGAT"
            "CGTCAGCCGTCAGCAACGAGTATTGTTGCTGACCCTTAACA"
            "ATCCCGCCGCACGTAATGCGCTAACTAATGCCCTGCTG", 5, GetParam()),
        std::set<std::string>({"AACAA"})
    );
}

TEST_P(TestFrequentWords, OutputAllMostFrequentKmer) {
    // This test dataset checks if your code correctly handles ties
    //  (i.e. your code actually outputs ALL “most frequent” kmers,
    //  and not just a single “most frequent” kmer).
    EXPECT_EQ(
        FrequentWords(
            "CCAGCGGGGGTTGATGCTCTGGGGGTCACAAGATTGCATTTTT"
            "ATGGGGTTGCAAAAATGTTTTTTACGGCAGATTCATTTAAA"
            "ATGCCCACTGGCTGGAGACATAGCCCGGATGCGCGTCTTTT"
            "ACAACGTATTGCGGGGTAAAATCGTAGATGTTTTAAAATAG"
            "GCGTAAC", 5, GetParam()),
        std::set<std::string>({"AAAAT", "GGGGT", "TTTTA"})
    );
}

TEST_P(TestFrequentWords, HandleEmptyPattern) {
    EXPECT_EQ(
        FrequentWords("", 3, GetParam()),
        std::set<std::string>({})
    );
    
    EXPECT_EQ(
        FrequentWords("ATCGTAGTCGCTAG", 0, GetParam()),
        std::set<std::string>({})
    );
}

TEST_P(TestFrequentWords, HandleOutOfBounds) {
    // k-mer can not be negative.    
    EXPECT_EQ(
        FrequentWords("ATCGTAGTCGCTAG", -4, GetParam()),
        std::set<std::string>({})
    );

    // k-mer equals to sequence length
    EXPECT_EQ(
        FrequentWords("ATG", 3, GetParam()),
        std::set<std::string>({"ATG"})
    );

    // k-mer is larger than sequence length    
    EXPECT_EQ(
        FrequentWords("ATG", 4, GetParam()),
        std::set<std::string>({})
    );
}


INSTANTIATE_TEST_SUITE_P(
    TestAllFrequentWords, TestFrequentWords,
    Values(AlgorithmEfficiency::Slow, AlgorithmEfficiency::Fast),
    [](const testing::TestParamInfo<TestFrequentWords::ParamType>& info) {
        switch (info.param)
        {
        case AlgorithmEfficiency::Slow:
            return "Slow";
        case AlgorithmEfficiency::Fast:
            return "Fast";
        default:
            return "Unknown";
        }
    }
);


TEST(TestPatternIndex, HandleNormalInput) {
    EXPECT_EQ(
        PatternIndex("ATGATGAAAATG", "ATG"),
        std::vector<size_t>({0, 3, 9})
    );
    
    EXPECT_EQ(
        PatternIndex("ATATAAAAGGGGATATA", "ATA"),
        std::vector<size_t>({0, 2, 12, 14})
    );
    
    EXPECT_EQ(
        PatternIndex("GATATATGCATATACTT", "ATAT"),
        std::vector<size_t>({1, 3, 9})
    );
}


TEST(TestPatternIndex, LongSequence) {
    EXPECT_EQ(
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
            "TCGACACCAGCCAGAGACACCGGACACCGACACCCCGAACACCAACACACCCGA", "ACACCA"),
        std::vector<size_t>({
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
        })
    );
}

TEST(TestPatternIndex, DonotReverseComplements) {
    // This dataset checks if your code is written correctly but is also taking into account
    // reverse complements, which we are not yet doing. Even though the reverse complement of
    // “ACAC” (which is “GTGT”) occurs in Genome, we only want to count occurrences of “ACAC”
    // specifically, which only occurs at index 4.
    EXPECT_EQ(
        PatternIndex("TTTTACACTTTTTTGTGTAAAAA", "ACAC"),
        std::vector<size_t>({4})
    );
}

TEST(TestPatternIndex, OffByOneErrorsOfStart) {
    // This dataset checks for off-by-one errors at the beginning of Genome. Notice that “AAA”
    // occurs at the very beginning of Genome, so if you were to miss the first kmer of Genome, your
    // code would output the following: 46 51 74
    EXPECT_EQ(
        PatternIndex(
            "AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAAT"
            "TAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAG"
            "CGCACAGACAGATAATAATTACAGAGTACACAACATCCAT", "AAA"),
        std::vector<size_t>({0, 46, 51, 74})
    );
}

TEST(TestPatternIndex, OffByOneErrorsOfEnd) {
    // This dataset checks for off-by-one errors at the endof Genome. Notice that “TTT” occurs
    // at the very end of Genome, so if you were to miss the last kmer of Genome, your code would
    // output the following: 88 93 98
    EXPECT_EQ(
        PatternIndex(
            "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTC"
            "ACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTT"
            "TATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", "TTT"),
        std::vector<size_t>({88, 92, 98, 132})
    );
}

TEST(TestPatternIndex, OverlapPattern) {
    // This test dataset checks if your code correctly handles cases where instances of Pattern
    // overlap in Genome. In this case, if you did not count overlaps, you would only find the first and
    // last instances of ATA (ATATATA and ATATATA). However, there is indeed a third
    // occurrence, where the other two overlap (ATATATA).
    EXPECT_EQ(
        PatternIndex("ATATATA", "ATA"),
        std::vector<size_t>({0, 2, 4})
    );
}

TEST(TestPatternIndex, EmptyInput) {
    EXPECT_EQ(
        PatternIndex("", "ATA"),
        std::vector<size_t>()
    );
    
    EXPECT_EQ(
        PatternIndex("ATATATA", ""),
        std::vector<size_t>()
    );
    
    EXPECT_EQ(
        PatternIndex("", ""),
        std::vector<size_t>()
    );

    EXPECT_EQ(
        PatternIndex("ATG", "ATGAGTA"),
        std::vector<size_t>()
    );
}

TEST(TestPatternIndex, NonNucleotide) {
    EXPECT_EQ(
        PatternIndex("ABCDEFG", "ABC"),
        std::vector<size_t>({0})
    );
}

TEST(TestFrequencyTable, NormalInput) {
    EXPECT_EQ(
        FrequencyTable("ACGTTTCACGTTTTACGG", 3),
        StrNumDict({
            {"ACG", 3},
            {"CGT", 2},
            {"GTT", 2},
            {"TTT", 3},
            {"TTC", 1},
            {"TCA", 1},
            {"CAC", 1},
            {"TTA", 1},
            {"TAC", 1},
            {"CGG", 1}
        })
    );

    EXPECT_EQ(
        FrequencyTable("ATGCAGTAAT", 3),
        StrNumDict({
            {"ATG", 1},
            {"TGC", 1},
            {"GCA", 1},
            {"CAG", 1},
            {"AGT", 1},
            {"GTA", 1},
            {"TAA", 1},
            {"AAT", 1}
        })
    );
}

TEST(TestFrequencyTable, EmptyInput) {
    EXPECT_EQ(
        FrequencyTable("ACGTTTCACGTTTTACGG", 0),
        StrNumDict()
    );

    EXPECT_EQ(
        FrequencyTable("", 3),
        StrNumDict()
    );

    EXPECT_EQ(
        FrequencyTable("", 0),
        StrNumDict()
    );

    EXPECT_EQ(
        FrequencyTable("ATGA", 6),
        StrNumDict()
    );

    EXPECT_EQ(
        FrequencyTable("ATG", 3),
        StrNumDict({
            {"ATG", 1}
        })
    );
}

class TestPatternToNumber : public TestWithParam<AlgorithmEfficiency> {
  // You can implement all the usual fixture class members here.
  // To access the test parameter, call GetParam() from class
  // TestWithParam<T>.
};

TEST_P(TestPatternToNumber, NormalInput) {
    EXPECT_EQ(
        PatternToNumber("A"),
        0b00
    );

    EXPECT_EQ(
        PatternToNumber("C"),
        0b01
    );

    EXPECT_EQ(
        PatternToNumber("G"),
        0b10
    );

    EXPECT_EQ(
        PatternToNumber("T"),
        0b11
    );

    EXPECT_EQ(
        PatternToNumber("ACGT"),
        0b00011011
    );

    EXPECT_EQ(
        PatternToNumber("ATG"),
        0b001110
    );

    EXPECT_EQ(
        PatternToNumber("TAA"),
        0b110000
    );

    EXPECT_EQ(
        PatternToNumber("AGT"),
        11
    );

    EXPECT_EQ(
        PatternToNumber("CTTCTCACGTACAACAAAATC"),
        2161555804173
    );
}

INSTANTIATE_TEST_SUITE_P(
    TestAllPatternToNumber, TestPatternToNumber,
    Values(AlgorithmEfficiency::Slow, AlgorithmEfficiency::Fast),
    [](const testing::TestParamInfo<TestFrequentWords::ParamType>& info) {
        switch (info.param)
        {
        case AlgorithmEfficiency::Slow:
            return "Slow";
        case AlgorithmEfficiency::Fast:
            return "Fast";
        default:
            return "Unknown";
        }
    }
);

TEST(TestNumberToPattern, NormalInput) {
    EXPECT_EQ(
        NumberToPatternBitwise(0b00, 1),
        "A"
    );
    EXPECT_EQ(
        NumberToPatternBitwise(0b01, 1),
        "C"
    );
    EXPECT_EQ(
        NumberToPatternBitwise(0b10, 1),
        "G"
    );
    EXPECT_EQ(
        NumberToPatternBitwise(0b11, 1),
        "T"
    );
    EXPECT_EQ(
        NumberToPatternBitwise(0b00011011, 4),
        "ACGT"
    );
    EXPECT_EQ(
        NumberToPatternBitwise(0b001110, 3),
        "ATG"
    );
    EXPECT_EQ(
        NumberToPatternBitwise(0b110000, 3),
        "TAA"
    );

    EXPECT_EQ(
        NumberToPatternBitwise(45, 4),
        "AGTC"
    );

    EXPECT_EQ(
        NumberToPatternBitwise(5353, 7),
        "CCATGGC"
    );

    EXPECT_EQ(
        NumberToPatternBitwise(7633, 10),
        "AAACTCTCAC"
    );
}

TEST(TestMaxMap, NormalInput) {
    EXPECT_EQ(
        MaxMap(
            StrNumDict({
                {"ATCT", 5},
                {"AGAAC", 6},
                {"GAT", 1}
            })
        ),
        6
    );

    EXPECT_EQ(
        MaxMap(
            StrNumDict({
                {"ATCT", 5},
                {"AGAA", 6},
                {"GATA", 1},
                {"GGGATA", 13},
                {"GAATTTA", 100},
                {"CATA", 61},
                {"GATAA", 14},
                {"TTATA", 1},
                {"GAAATA", 13},
                {"GCAATA", 7},
                {"GAGATA", 15},
                {"GACGTTA", 100},
            })
        ),
        100
    );
}

TEST(TestMaxMap, EmptyInput) {
    EXPECT_THROW(
        MaxMap(StrNumDict()),
        std::runtime_error
    );
}

} // namespace