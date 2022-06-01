#include <string>
#include <set>
#include <vector>
#include <iostream>

#include "gtest/gtest.h"

#include "pattern.h"
#include "utils.h"

namespace {

using namespace bioutils::algorithms;
using namespace bioutils::utils;
using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::ValuesIn;

typedef std::unordered_map<std::string, uint> StrNumDict;

class TestPatternCount: public TestWithParam<AlgorithmEfficiency> {};

TEST_P(TestPatternCount, NormalInput) {
    EXPECT_EQ(
        PatternCount("GCGCG", "GCG", GetParam()),
        2
    );

    /*
        This dataset just checks if you’re correctly counting. It is the “easiest” test. Notice that all                               
        occurrences of CG in Text (A​CG​TA​CG​TA​CG​T) are away from the very edges (so your code                             
        won’t fail on off­by­one errors at the beginning or at the end of Text) and that none of the                                     
        occurrences of Pattern overlap (so your code won’t fail if you fail to account for overlaps).
    */
    EXPECT_EQ(
        PatternCount("ACGTACGTACGT", "CG", GetParam()),
        3
    );

    /*
        This dataset checks if your code correctly handles cases where there is an occurrence of                             
        Pattern at the very beginning of Text. Note that there are no overlapping occurrences of Pattern                               
        (i.e. AAAA), and there is no occurrence of Pattern at the very end of Text, so assuming your                                   
        code passed Test Dataset 1, this test would only check for off­by­one errors at the beginning of                                 
        Text.
    */
    EXPECT_EQ(
        PatternCount(
            "AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCT"
            "GCCGTGAGTAAATTAAATTTTATTGACTTAGGTCACT"
            "AAATACTTTAACCAATATAGGCATAGCGCACAGACAG"
            "ATAATAATTACAGAGTACACAACATCCAT", "AAA", GetParam()),
        4
    );

    /*
        This dataset checks if your code correctly handles cases where there is an occurrence of                             
        Pattern at the very end of Text. Note that there are no overlapping occurrences of Pattern (i.e.                                 
        AAAA), and there is no occurrence of Pattern at the very beginning of Text, so assuming your                                 
        code passed Test Dataset 2, this test would only check for off­by­one errors at the end of Text.
    */
    EXPECT_EQ(
        PatternCount(
            "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGG"
            "CCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAA"
            "GCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCG"
            "TTAACCACGTTCTGTGCCGACTTT", "TTT", GetParam()),
        4
    );

    /*
        This test dataset checks if your code is also counting occurrences of the Reverse                           
        Complement of Pattern (which would have an output of 4), which is out of the scope of this                                   
        problem (that will come up later in the chapter). Your code should only be looking for perfect                                 
        matches of Pattern in Text at this point.
    */
    EXPECT_EQ(
        PatternCount("GGACTTACTGACGTACG", "ACT", GetParam()),
        2
    );

    /*
        This dataset checks if your code correctly handles cases where occurrences of Pattern                         
        overlap. For example, any occurrence of the string “CCC” should count as 2 occurrences of                             
        “CC” (​CC​C and C​CC​). In this dataset, there are 5 occurrences of CC including overlaps                             
        (AT​CC​GAT​CCC​ATG​CCC​ATG).
    */
    EXPECT_EQ(
        PatternCount("ATCCGATCCCATGCCCATG", "CC", GetParam()),
        5
    );

    // This is the final test that we run your code on: the Full Dataset.
    EXPECT_EQ(
        PatternCount(
            "CTGTTTTTGATCCATGATATGTT"
            "ATCTCTCCGTCATCAGAAGAACA"
            "GTGACGGATCGCCCTCTCTCTTG"
            "GTCAGGCGACCGTTTGCCATAAT"
            "GCCCATGCTTTCCAGCCAGCTCT"
            "CAAACTCCGGTGACTCGCGCAGG"
            "TTGAGTA", "CTC", GetParam()),
        9
    );
}

INSTANTIATE_TEST_SUITE_P(
    TestAllPatternCount, TestPatternCount,
    Values(AlgorithmEfficiency::Slow, AlgorithmEfficiency::Fast, AlgorithmEfficiency::Faster),
    [](const testing::TestParamInfo<TestPatternCount::ParamType>& info) {
        switch (info.param)
        {
        case AlgorithmEfficiency::Slow:
            return "Slow";
        case AlgorithmEfficiency::Fast:
            return "Fast";
        case AlgorithmEfficiency::Faster:
            return "Faster";
        default:
            return "Unknown";
        }
    }
);

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
    Values(
        AlgorithmEfficiency::Slow,
        AlgorithmEfficiency::Fast,
        AlgorithmEfficiency::Faster,
        AlgorithmEfficiency::Fastest),
    [](const testing::TestParamInfo<TestFrequentWords::ParamType>& info) {
        switch (info.param)
        {
        case AlgorithmEfficiency::Slow:
            return "FrequentWordsSlow";
        case AlgorithmEfficiency::Fast:
            return "FrequentWordsByPerfectHash";
        case AlgorithmEfficiency::Faster:
            return "FrequentWordsBySorting";
        case AlgorithmEfficiency::Fastest:
            return "FrequentWordsByStdHash";
        default:
            return "Unknown";
        }
    }
);

typedef std::set<std::string> (*FrequentWordsWithMismatchesFuncPtr)(std::string_view, int, int, bool rev);
class TestFrequentWordsWithMismatches : public TestWithParam<FrequentWordsWithMismatchesFuncPtr> {};

TEST_P(TestFrequentWordsWithMismatches, HandleNormalInput) {
    FrequentWordsWithMismatchesFuncPtr fun = GetParam();

    EXPECT_EQ(
        fun("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1, false),
        std::set<std::string>({"GATG", "ATGC", "ATGT"})
    );

    /*
        Text contains partial and complete matches for the most frequent word.
    */
    EXPECT_EQ(
        fun("AGGT", 2, 1, false),
        std::set<std::string>({"GG"})
    );

    EXPECT_EQ(
        fun("AGGGT", 2, 0, false),
        std::set<std::string>({"GG"})
    );

    /*
        Text has multiple most frequent words
    */
    EXPECT_EQ(
        fun("AGGCGG", 3, 0, false),
        std::set<std::string>({"AGG", "GGC", "GCG", "CGG"})
    );
}

TEST_P(TestFrequentWordsWithMismatches, ReverseComplements) {
    FrequentWordsWithMismatchesFuncPtr fun = GetParam();

    EXPECT_EQ(
        fun("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1, true),
        std::set<std::string>({"ATGT", "ACAT"})
    );

    /*
        Text contains partial and completes matches for the most frequent word.
    */
    EXPECT_EQ(
        fun("AAAAAAAAAA", 2, 1, true),
        std::set<std::string>({"AT", "TA"})
    );

    /*
        This dataset makes sure that your code is not accidentally swapping k and d
    */
    EXPECT_EQ(
        fun("AGTCAGTC", 4, 2, true),
        std::set<std::string>({"AATT", "GGCC"})
    );

    /*
        This dataset makes sure you are finding k-mers in both Text and
        the reverse complement of Text.
    */
    EXPECT_EQ(
        fun("AATTAATTGGTAGGTAGGTA", 4, 0, true),
        std::set<std::string>({"AATT"})
    );

    /*
        This dataset first checks that k-mers with exactly d mismatches are being found. Then,
        it checks that k-mers with less than d mismatches are being allowed (i.e. you are not only allowing
        k-mers with exactly d mismatches). Next, it checks that you are not returning too few k-mers. Last,
        it checks that you are not returning too many k-mers.
    */
    EXPECT_EQ(
        fun("ATA", 3, 1, true),
        std::set<std::string>({
            "AAA", "AAT", "ACA", "AGA", "ATA",
            "ATC", "ATG", "ATT", "CAT", "CTA",
            "GAT", "GTA", "TAA", "TAC", "TAG",
            "TAT", "TCT", "TGT", "TTA", "TTT"})
    );

    /*
        This dataset checks that your code is looking at both Text and its reverse complement
        (i.e. not just looking at Text, and not just looking at the reverse complement of Text, but looking at
        both).
    */
    EXPECT_EQ(
        fun("AAT", 3, 0, true),
        std::set<std::string>({"AAT", "ATT"})
    );

    /*
        This dataset checks that your code correctly delimiting your output (i.e. using spaces)
        and verifies that your k-mers are actually of length k.
    */
    EXPECT_EQ(
        fun("TAGCG", 2, 1, true),
        std::set<std::string>({"CA", "CC", "GG", "TG"})
    );
}

INSTANTIATE_TEST_SUITE_P(
    TestAllFrequentWordsWithMismatches, TestFrequentWordsWithMismatches,
    Values(
        FrequentWordsWithMismatches,
        FrequentWordsWithMismatchesBySorting),
    [](const ::testing::TestParamInfo<TestFrequentWordsWithMismatches::ParamType>& info) {
        if (info.param == FrequentWordsWithMismatches) {
            return "FrequentWordsWithMismatches";
        } else if (info.param == FrequentWordsWithMismatchesBySorting) {
            return "FrequentWordsWithMismatchesBySorting";
        } else {
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

TEST(TestFrequencyArray, NormalInput) {
    EXPECT_EQ(
        FrequencyArray("ACGCGGCTCTGAAA", 2),
        std::vector<uint>(
            {2, 1, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 1, 1, 0})
    );

    /*
        This dataset checks if you have an off-­by-­one error at the end of Text
        (i.e. you are not counting the last kmer in Text). There are three
        instances of AA (​AA​AAC, A​AA​AC, and AA​AA​C),  but  there  is  one
        instance  of  AC  at  the  end  (AAA​AC​).
    */
    EXPECT_EQ(
        FrequencyArray("AAAAC", 2),
        std::vector<uint>(
            {3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
        )
    );

    /*
        This dataset checks if you have an off-­by-­one error at the beginning
        of Text (i.e. you are not counting the first kmer in Text). There are
        two instances of AA (TT​AA​A and TTA​AA​), but there  is  one  instance
        of  TTA  (​TTA​AA)  and  one  instance  of  TAA  (T​TAA​A).
    */
    EXPECT_EQ(
        FrequencyArray("TTAAA", 2),
        std::vector<uint>(
            {2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1}
        )
    );

    /*
        This dataset checks if your code actually increments each count, or if your
        code instead just sets the count equal to one each time. In other words,
        this dataset checks if your code is doing  something like ​array[kmer] = 1​
        instead  of  ​array[kmer] += 1​. 
    */
    EXPECT_EQ(
        FrequencyArray("AAA", 2),
        std::vector<uint>(
            {2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
        )
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

    // Max hashable length
    EXPECT_EQ(
        PatternToNumber("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
        0b1111111111111111111111111111111111111111111111111111111111111111
    );

    EXPECT_THROW(
        PatternToNumber("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTC"),
        std::runtime_error
    );
}

INSTANTIATE_TEST_SUITE_P(
    TestAllPatternToNumber, TestPatternToNumber,
    Values(AlgorithmEfficiency::Slow, AlgorithmEfficiency::Fast),
    [](const testing::TestParamInfo<TestPatternToNumber::ParamType>& info) {
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

TEST(TestMaxArray, NormalInput) {
    EXPECT_EQ(
        MaxArray(std::vector<uint>{1, 2, 3, 4, 5}),
        5
    );

    EXPECT_EQ(
        MaxArray(std::vector<uint>{2, 5, 6, 3, 2, 5, 12, 9, 10}),
        12
    );
}

TEST(TestMaxArray, EmptyInput) {
    EXPECT_THROW(
        MaxArray(std::vector<uint>{}),
        std::runtime_error
    );
}

class TestFindClumps: public TestWithParam<AlgorithmEfficiency> {};

TEST_P(TestFindClumps, NormalInput) {
    EXPECT_EQ(
        FindClumps(
            "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAA"
            "GAGAAGAGGAAACACGACACGACATTGCGACATAATGTAC"
            "GAATGTAATGTGCCTATGGC", 5, 75, 4, GetParam()),
        std::set<std::string>({"CGACA", "GAAGA", "AATGT"})
    );

#if 0
    EXPECT_EQ(
        FindClumps(
            "GCGGTTATGCACCGTTCAAATTAGCAAACCACTAAGCGAC"
            "GTAGTCTGGATTGATTTCTCCCTACCAGTGACCCAAGACG"
            "CGTTAGTGAGTTAAGTTCATATCCAGTACCTGCCGCCCTC"
            "TGTACTTGGGCGTCCGATTCGCATGCTTACTCAGGTGGAG"
            "GACACGATAATCTGATTAAACTGAGCTAAACCAGGTGGAA"
            "CCAGAAACCAGGTGGGGAGTCTCGCTTCAAGCCGTTCTTG"
            "CGATCAAACCAGGTGGTCCATTATGAAACCAGGTGGCTAA"
            "ACCAGGTGGTCCAGATCCTCGAATGATGTCGGTGCACATC"
            "AAAACCAGGTGGGGTGGTGGAACGTAAAACCAGGTGGCAT"
            "AAACCAGGTGGGCCGGTTCGTAAACCAGGTGAAACCAGGT"
            "GGGGTGGAAACCAGGTGGGTTACAAATTACGTTGAGATGG"
            "CCCAAACCAGGTGGTGGGCTTCACCCATGTCAACAAACCA"
            "CCCTATGGAACTAAACCAGGTGGAACCAGGTGGTGAAGGC"
            "TTATCCTCAGGAAAAACCAGGTGGAGGTGGTGAAATAAAA"
            "CCAGGTGGACCAGGTGGATAACCCTCGCCTCGCTTCTCAA"
            "CCGAGACCTGGATAAACCAGGTGGGGTGGTCCACCGATTT"
            "TTGAGACACTAGAAACCAGGTGGGCGGGGAAACCAGGTGG"
            "CAAACCAGGTGGGGTGGACGGAAACCAGGTGGATATGTCA"
            "TAAAACCAAACCAGGTGGTGCACCCCCATGGTGTGTCTTA"
            "TCCGTGCGTATAAACCAGGTGGTCGCACGGCTTCCACTTG"
            "CTGAGAATAGGCCCGCAGGGTCAGTGCCATGCCCTCCGTC"
            "ACTCGATATGTGTTGTAAGAGTGGTTACCCCTTCATTGAA"
            "GTCGCCCACAGCCCCACCTGCATTGCTAGACTATCACCCT"
            "ACAGTAGGCCTTTTCGCCTTCTTCAAGCAGCAATCTCTTA"
            "TCCGCGGATGGGCGCGGCGAGCGTGGCGTCCCCGAACATT"
            "TTTACCTAACGTGTTTTGTTGGCCGCAAGCCTTCCCTCTA"
            "GTCCACCTCAGCCATTCAGCCTAGTAGCTTTCAAGCCGAG"
            "CCTTCCATATCTAATGGACCGTCCAGAATTTCACACGTTT"
            "CACAGGGCTGTGTTCGACCGCCCGTAATGCTGTTTCACAG"
            "GCGATCGCCTTGCGGTTTTTTCACAGATCGCAGCCGATGG"
            "ACATGCCAACTCGATTTTCACAGAGTTTTTCACAGCGGTT"
            "TCACAGCACAGCAGTGATTGTTTCACAGCAATTTTCACTT"
            "TCACAGGGGCCCTTTTCACAGCTCAGGGCTCTTTTCACTT"
            "TCACAGTTTCACAGCGCTCCTTTCACAGAGCGGGGAAATT"
            "TAAGGGAACACTCAAGGGAACAAGGGAACACACAAAGGGA"
            "ACACAACACAACACATAAGGGAACACTTTCACAGAACACA"
            "AAAGTCCGAAATCATCAGCGGCGAAGGGATTTCACAGACA"
            "GACACTTTCACAGCGCATTTCACAGATACGTACTTTCACA"
            "GGCGTACTTTCACAGACTTTCACAGAGGACAAGCTCAATT"
            "TTCACAGACAGGCTGGATAAATTTCACAGCGGTAAGGGTT"
            "TCACAGCACACATAAGGGAACACGAATTTCACAGCAGGGA"
            "ACACCTCTACGAGTAATCTATTACTCTACCTACTGAAGGG"
            "AACACACCGAAGACCTACTATTACCTATTACTCTTAAAGG"
            "GAACACATTACAAGGGAACACACTCTCTCGTCATATCTCA"
            "CCTCTCTATTACTCTTAAGGGAACACCTTCTCGATCAACC"
            "TATTACTCTATGGAGATAGAGATATTCCAGACATATGGAG"
            "ATAACATGGAGATATGGAGATAATGGAGATGGAGATAGCT"
            "CTTATATTTATCCTATGGAGATATGATACTATTAATGGAG"
            "ATAATTCTAATGGAGATATAATTACTCTAAGAGGATGGGA"
            "TCTCGGGCTATTACTCTAATGGAGATAAGCACTATTACTC"
            "TAGGAAATGGAGATATGTCAATGGAGATATGTAATGGAGA"
            "TAGAGGGAGATGGAGTCGCCATTTCATAATCGCCATTTCA"
            "TAGTTCAGGAATCGCCATTTCCGCCATTTCTAAGATGGAG"
            "TCGCCATTTCTACGTATGGAGATAGGATCGCCATTTCATA"
            "CGACCCGTTGGATATCGCCATTTCCTCGCCATTTCTGGTG"
            "ACATTTCTCGCCATTTCATTTCTGGAGATAGATGGATCTC"
            "GCCATTTCATAGGAATCGCCATTTCCACGTAGGGGGGGCC"
            "ACAATCCGTAGGTCGGAATTCAGACTCGCCATTTCCCATC"
            "GCCATTTCTTCACCTGTATGCCGATCCCTTCGCCATTTCT"
            "CATGGAGATAACTCTCTCTCGCCATTTCTCGCCATTTCCA"
            "TTTCACTCTCATTCGCCATCGCCATTTCCATTCGCCATTT"
            "CATCGCCATTTCTTCAGGATAAGATATCGCCATTTCGACT"
            "CTCATTCGCATACTGACTCTCATTCTCATCTCGCCATTTC"
            "TCATCTGACTCTCATCCTGGGGGAAACTTGCGACTCTCAT"
            "CACACTTCCGTCGACTCTCATACTGGCGGATAGCATAGGA"
            "GCCATTTAAAGACTCTCATTCTCATTCGAGACTCTCATTC"
            "AAATCCTACGAGGACTCTCATATAGACTCTCATATCATTA"
            "CGAGGACTCTCATATACGAGCCATGCATGTGGCGACGACT"
            "CTCATCTACGAGCCATGCAAGCAGAATCTACGAGCGACTC"
            "TCATTACGAGCCATGTGACCGTACGAGCCATGCATGCATG"
            "CCATGCTGACTCTCATCGAGTACGAGCCATGGAAGTTCTT"
            "GTTGGTTCGTAGCCCAAGAGCTGAAGTTACGAGCCTACGA"
            "GCCATGAAGTTACTTTTACGAGCCATGAAGCTTACGATAC"
            "GAGCCATGCGAGCCATGCATCCGCGCTACGAGCCATGTTC"
            "CAGTACGAGCCATGTTAGTTGCTGAAGTTAAGTTTGGCGC"
            "TGAAGTTTGTACGAGCCATGTGCCCGCTGAAGTTTGTTGT"
            "ACGAGCCATGCATGCTGAAGTTAATGGCTGAAGTTAGCGT"
            "TTGCGGGCAGATCCTCATTCTACGATACGAGCCATGCCAT"
            "GCAGCTGAAGTTAAGTTGGGTTACGAGCCATGCGAGCCAT"
            "GTGAAGTACGAGCCATGCTGGCTGAAGTTGTTTGTGCTGC"
            "TGAAGTTGCTCTTGTCTCTAGCTGAAGTTGCCAACAGGGC"
            "TGAAGCTGAAGTTTAAGCTGAAGTTGCGAGCAGGCTGAAG"
            "TTATCGGATTGGGGCTGAAGTTCAACCTCCCGTCCCCCCA"
            "CACTATATTCCCGTCCCCCCCCGCGCACGCGCCGTCTCCC"
            "GTCCCCCCTATCCCGTGCGCACGCGACGCGATCCCGTCCC"
            "CCCAGAGTGCGCGCACGCGTCCCCCTTCCCGTCCCCCTCT"
            "CCCGGGCGCACGCGTCGCTCAACATTTCCGCGCACGCGTC"
            "GCGCACGCGGGCGCACGCGGGTCCCGTCCCCCCCCCTCTT"
            "CGGCGCACGCGGAATTCCCGTCGCGCACGCGTCCCGTCCC"
            "GCGCACGCGTCGCGCACGCGACTGCCCTAACCAACAGTGC"
            "GCACGCGCCGGTAACCCGGTAACCCGGTAACCGCGCACGC"
            "GGGCGCACGCGCGTAACCCGCGCACGCGCCGCGCACGCGG"
            "CCCGGTTCCCGTCCCCCCCGGTAACCCGGTAACTCCCGTC"
            "CCCCGTAACCCGGTGCGCACGCGCCCGGCGCACGCGGAGC"
            "GCACGCGCCCCCCCCGGTAATAGCGCACGCGCCCGGGCGC"
            "ACGCGCCCGGTAACCCGGTAACCCGGGCGCGCGCACGCGG"
            "CGGCGCACGCGGCGCACGCGGCGCACGCG", 11, 566, 18, GetParam()),
        std::set<std::string>({"AAACCAGGTGG"})
    );
#endif

    /*
        This dataset makes sure that your code only counts kmers that fall COMPLETELY
        within a given L-­window. For example, take the 4-­window starting at index 4
        (AAAA​CGTC​GAAAAA). One might think that the 2-­mer “CG” occurs twice in this window
        since the first letter of the second occurrence happens at the very end of the window. However,
        since the second occurrence of “CG” does not fall entirely in this 4-­window, it does not count.
        Thus,  the  only  result  is  “AA”.   
    */
    EXPECT_EQ(
        FindClumps("AAAACGTCGAAAAA", 2, 4, 2, GetParam()),
        std::set<std::string>({"AA"})
    );

    /*
        This dataset checks if your code has an off-­by-­one error when checking kmers within an
        L-­window. Notice that, for each 1-­mer (A, C, G, and T), there are 3 nucleotides between the first
        and second occurrence. In other words, each nucleotide occurs twice in a specific 5-­window:
        once at the beginning of the 5-­window, and once at the end​.
    */
    EXPECT_EQ(
        FindClumps("ACGTACGT", 1, 5, 2, GetParam()),
        std::set<std::string>({"A", "C", "G", "T"})
    );

    /*
        This dataset checks if your code is correctly handling overlapping kmers.
    */
    EXPECT_EQ(
        FindClumps(
            "CCACGCGGTGTACGCTGCAAAAAGCCTTGCTGAATCAAAT"
            "AAGGTTCCAGCACATCCTCAATGGTTTCACGTTCTTCGCC"
            "AATGGCTGCCGCCAGGTTATCCAGACCTACAGGTCCACCA"
            "AAGAACTTATCGATTACCGCCAGCAACAATTTGCGGTCCA"
            "TATAATCGAAACCTTCAGCATCGACATTCAACATATCCAG"
            "CG", 3, 25, 3, GetParam()),
        std::set<std::string>({"AAA", "CAG", "CAT", "CCA", "GCC", "TTC"})
    );
}

INSTANTIATE_TEST_SUITE_P(
    TestAllFindClumps, TestFindClumps,
    Values(
        AlgorithmEfficiency::Slow,
        AlgorithmEfficiency::Fast,
        AlgorithmEfficiency::Faster,
        AlgorithmEfficiency::Fastest),
    [](const testing::TestParamInfo<TestFindClumps::ParamType>& info) {
        switch (info.param)
        {
        case AlgorithmEfficiency::Slow:
            return "FindClumpsRaw";
        case AlgorithmEfficiency::Fast:
            return "FindClumpsRaw2";
        case AlgorithmEfficiency::Faster:
            return "FindClumpsBetterWithPerfectHash";
        case AlgorithmEfficiency::Fastest:
            return "FindClumpsBetterWithStdHash";
        default:
            return "Unknown";
        }
    }
);

TEST(TestFindMinimumSkew, NormalInput) {
    EXPECT_EQ(
        FindMinimumSkew("CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG"),
        std::vector<size_t>({53, 97})
    );

    EXPECT_EQ(
        FindMinimumSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"),
        std::vector<size_t>({11, 24})
    );

    /*
        This dataset checks if your code’s indexing is off. Specifically, it verifies that your code
        is not returning an index 1 too high (i.e. 4) or 1 too low (i.e. 2). 
    */
    EXPECT_EQ(
        FindMinimumSkew("ACCG"),
        std::vector<size_t>({3})
    );

    /*
        This dataset checks to see if your code is missing the last symbol of ​Genome.
    */
    EXPECT_EQ(
        FindMinimumSkew("ACCC"),
        std::vector<size_t>({4})
    );

    /*
        This dataset makes sure you’re not accidentally finding the maximum skew instead of the
        minimum skew.
    */
    EXPECT_EQ(
        FindMinimumSkew("CCGGGT"),
        std::vector<size_t>({2})
    );

    /*
        First, this dataset checks if you are only finding 1 index (and not multiple indices). Then,
        it checks if you are using a delimiter to separate your indices (ideally a space character).
    */
    EXPECT_EQ(
        FindMinimumSkew("CCGGCCGG"),
        std::vector<size_t>({2, 6})
    );
}

TEST(TestHammingDistance, NormalInput) {
    EXPECT_EQ(
        HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC"),
        3
    );

    /*
        This dataset checks if your code isn’t keeping count (i.e. returns ‘0’ when the answer is
        clearly nonzero) or if your code returns a negative value, which is impossible.
    */
    EXPECT_EQ(
        HammingDistance("AAAA", "TTTT"),
        4
    );

    /*
        This dataset checks if your code is finding Edit Distance (which would be 2) instead of
        Hamming Distance.
    */
    EXPECT_EQ(
        HammingDistance("ACGTACGT", "TACGTACG"),
        8
    );

    /*
        This dataset checks if your code is returning the number of matches (2) instead of the
        number of mismatches (6).
    */
    EXPECT_EQ(
        HammingDistance("ACGTACGT", "CCCCCCCC"),
        6
    );

    /*
        This dataset checks if your code works on a dataset where the two input strings have no
        matches.
    */
    EXPECT_EQ(
        HammingDistance("ACGTACGT", "TGCATGCA"),
        8
    );

    /*
        This dataset checks if you have an off­by­one error at the beginning (i.e. you are starting
        at the second character of the strings instead of the first character).
    */
    EXPECT_EQ(
        HammingDistance(
            "GATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACT",
            "AATAGCAGCTTCTCAACTGGTTACCTCGTATGAGTAAATTAGGTCATTATTGACTCAGGTCACTAACGTCT"),
        15
    );

    /*
        This dataset checks if you have an off­by­one error at the end (i.e. you are ending at the
        second­to­last character of the strings instead of the last character).
    */
    EXPECT_EQ(
        HammingDistance(
            "AGAAACAGACCGCTATGTTCAACGATTTGTTTTATCTCGTCACCGGGATATTGCGGCCACTCATCGGTCAGTTGATTACGCAGGGCGTAAATCGCCAGAATCAGGCTG",
            "AGAAACCCACCGCTAAAAACAACGATTTGCGTAGTCAGGTCACCGGGATATTGCGGCCACTAAGGCCTTGGATGATTACGCAGAACGTATTGACCCAGAATCAGGCTC"),
        28
    );
}

TEST(TestPatternIndexApproximate, NormalInput) {
    EXPECT_EQ(
        PatternIndexApproximate(
            "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC",
            "ATTCTGGA", 3),
        std::vector<size_t>({6, 7, 26, 27, 78})
    );

    /*
        This dataset checks if you are only counting instances where the number of mismatches is
        exactly equal to d (i.e. ignoring instances where mismatch < d). 
    */
    EXPECT_EQ(
        PatternIndexApproximate(
            "TTTTTTAAATTTTAAATTTTTT", "AAA", 2),
        std::vector<size_t>({4, 5, 6, 7, 8, 11, 12, 13, 14, 15})
    );

    /*
        This dataset checks if your code has an off-by-one error at the beginning of Text (i.e. your 
        code is not checking the the leftmost substring of Text).
    */
    EXPECT_EQ(
        PatternIndexApproximate(
            "GAGCGCTGGGTTAACTCGCTACTTCCCGACGAGCGCTGTGGCGCAAATTGGCGATGA"
            "AACTGCAGAGAGAACTGGTCATCCAACTGAATTCTCCCCGCTATCGCATTTTGATGC"
            "GCGCCGCGTCGATT", "GAGCGCTGG", 2),
        std::vector<size_t>({0, 30, 66})
    );

    /*
        This dataset checks if your code has an off-by-one error at the end of Text (i.e. your code
        is not checking the the rightmost substring of Text).
    */
    EXPECT_EQ(
        PatternIndexApproximate(
            "CCAAATCCCCTCATGGCATGCATTCCCGCAGTATTTAATCCTTTCATTCTGCATATAA"
            "GTAGTGAAGGTATAGAAACCCGTTCAAGCCCGCAGCGGTAAAACCGAGAACCATGA"
            "TGAATGCACGGCGATTGCGCCATAATCCAAACA", "AATCCTTTCA", 3),
        std::vector<size_t>({3, 36, 74, 137})
    );

    /*
        This  dataset  checks  if  your  code  is  correctly  accounting  for  overlapping  instances  of
        Pattern in Text.
    */
    EXPECT_EQ(
        PatternIndexApproximate(
            "CCGTCATCCGTCATCCTCGCCACGTTGGCATGCATTCCGTCATCCCGTCAGGCATACT"
            "TCTGCATATAAGTACAAACATCCGTCATGTCAAAGGGAGCCCGCAGCGGTAAAACC"
            "GAGAACCATGATGAATGCACGGCGATTGC", "CCGTCATCC", 3),
        std::vector<size_t>({0, 7, 36, 44, 48, 72, 79, 112})
    );

    /*
        This  dataset  checks  if  you  are  only  counting  instances  of  Pattern  with  less  than  d
        mismatches (as opposed to instances of Pattern with less than or equal to d mismatches).
    */
    EXPECT_EQ(
        PatternIndexApproximate(
            "AAAAAA", "TTT", 3),
        std::vector<size_t>({0, 1, 2, 3})
    );

    /*
        This dataset checks if your code works with input where d = 0 (i.e. only perfect matches
        are allowed).
    */
    EXPECT_EQ(
        PatternIndexApproximate(
            "CCACCT", "CCA", 0),
        std::vector<size_t>({0})
    );
}

TEST(TestNeighbors, NormalInput) {

    auto test1 =  std::set<std::string>({
        "CCG",
        "TCG",
        "GCG",
        "AAG",
        "ATG",
        "AGG",
        "ACA",
        "ACC",
        "ACT",
        "ACG"
    });

    EXPECT_EQ(
        NeighborsRecursive("ACG", 1),
        test1
    );

    auto test2 = std::set<std::string>({
        "AACAATAT", "ACAAATAT", "ACCAAAAT", "ACCAACAT", "ACCAAGAT", "ACCAATAA", "ACCAATAC", "ACCAATAG",
        "ACCAATAT", "ACCAATCT", "ACCAATGT", "ACCAATTT", "ACCACTAT", "ACCAGTAT", "ACCATTAT", "ACCCATAT",
        "ACCGATAT", "ACCTATAT", "ACGAATAT", "ACTAATAT", "AGCAATAT", "ATCAATAT", "CACAATAT", "CCAAATAT",
        "CCCAAAAT", "CCCAACAT", "CCCAAGAT", "CCCAATAA", "CCCAATAC", "CCCAATAG", "CCCAATAT", "CCCAATCT",
        "CCCAATGT", "CCCAATTT", "CCCACTAT", "CCCAGTAT", "CCCATTAT", "CCCCATAT", "CCCGATAT", "CCCTATAT",
        "CCGAATAT", "CCTAATAT", "CGCAATAT", "CTCAATAT", "GACAATAT", "GCAAATAT", "GCCAAAAT", "GCCAACAT",
        "GCCAAGAT", "GCCAATAA", "GCCAATAC", "GCCAATAG", "GCCAATAT", "GCCAATCT", "GCCAATGT", "GCCAATTT",
        "GCCACTAT", "GCCAGTAT", "GCCATTAT", "GCCCATAT", "GCCGATAT", "GCCTATAT", "GCGAATAT", "GCTAATAT",
        "GGCAATAT", "GTCAATAT", "TAAAATAT", "TACAAAAT", "TACAACAT", "TACAAGAT", "TACAATAA", "TACAATAC",
        "TACAATAG", "TACAATAT", "TACAATCT", "TACAATGT", "TACAATTT", "TACACTAT", "TACAGTAT", "TACATTAT",
        "TACCATAT", "TACGATAT", "TACTATAT", "TAGAATAT", "TATAATAT", "TCAAAAAT", "TCAAACAT", "TCAAAGAT",
        "TCAAATAA", "TCAAATAC", "TCAAATAG", "TCAAATAT", "TCAAATCT", "TCAAATGT", "TCAAATTT", "TCAACTAT",
        "TCAAGTAT", "TCAATTAT", "TCACATAT", "TCAGATAT", "TCATATAT", "TCCAAAAA", "TCCAAAAC", "TCCAAAAG",
        "TCCAAAAT", "TCCAAACT", "TCCAAAGT", "TCCAAATT", "TCCAACAA", "TCCAACAC", "TCCAACAG", "TCCAACAT",
        "TCCAACCT", "TCCAACGT", "TCCAACTT", "TCCAAGAA", "TCCAAGAC", "TCCAAGAG", "TCCAAGAT", "TCCAAGCT",
        "TCCAAGGT", "TCCAAGTT", "TCCAATAA", "TCCAATAC", "TCCAATAG", "TCCAATAT", "TCCAATCA", "TCCAATCC",
        "TCCAATCG", "TCCAATCT", "TCCAATGA", "TCCAATGC", "TCCAATGG", "TCCAATGT", "TCCAATTA", "TCCAATTC",
        "TCCAATTG", "TCCAATTT", "TCCACAAT", "TCCACCAT", "TCCACGAT", "TCCACTAA", "TCCACTAC", "TCCACTAG",
        "TCCACTAT", "TCCACTCT", "TCCACTGT", "TCCACTTT", "TCCAGAAT", "TCCAGCAT", "TCCAGGAT", "TCCAGTAA",
        "TCCAGTAC", "TCCAGTAG", "TCCAGTAT", "TCCAGTCT", "TCCAGTGT", "TCCAGTTT", "TCCATAAT", "TCCATCAT",
        "TCCATGAT", "TCCATTAA", "TCCATTAC", "TCCATTAG", "TCCATTAT", "TCCATTCT", "TCCATTGT", "TCCATTTT",
        "TCCCAAAT", "TCCCACAT", "TCCCAGAT", "TCCCATAA", "TCCCATAC", "TCCCATAG", "TCCCATAT", "TCCCATCT",
        "TCCCATGT", "TCCCATTT", "TCCCCTAT", "TCCCGTAT", "TCCCTTAT", "TCCGAAAT", "TCCGACAT", "TCCGAGAT",
        "TCCGATAA", "TCCGATAC", "TCCGATAG", "TCCGATAT", "TCCGATCT", "TCCGATGT", "TCCGATTT", "TCCGCTAT",
        "TCCGGTAT", "TCCGTTAT", "TCCTAAAT", "TCCTACAT", "TCCTAGAT", "TCCTATAA", "TCCTATAC", "TCCTATAG",
        "TCCTATAT", "TCCTATCT", "TCCTATGT", "TCCTATTT", "TCCTCTAT", "TCCTGTAT", "TCCTTTAT", "TCGAAAAT",
        "TCGAACAT", "TCGAAGAT", "TCGAATAA", "TCGAATAC", "TCGAATAG", "TCGAATAT", "TCGAATCT", "TCGAATGT",
        "TCGAATTT", "TCGACTAT", "TCGAGTAT", "TCGATTAT", "TCGCATAT", "TCGGATAT", "TCGTATAT", "TCTAAAAT",
        "TCTAACAT", "TCTAAGAT", "TCTAATAA", "TCTAATAC", "TCTAATAG", "TCTAATAT", "TCTAATCT", "TCTAATGT",
        "TCTAATTT", "TCTACTAT", "TCTAGTAT", "TCTATTAT", "TCTCATAT", "TCTGATAT", "TCTTATAT", "TGAAATAT",
        "TGCAAAAT", "TGCAACAT", "TGCAAGAT", "TGCAATAA", "TGCAATAC", "TGCAATAG", "TGCAATAT", "TGCAATCT",
        "TGCAATGT", "TGCAATTT", "TGCACTAT", "TGCAGTAT", "TGCATTAT", "TGCCATAT", "TGCGATAT", "TGCTATAT",
        "TGGAATAT", "TGTAATAT", "TTAAATAT", "TTCAAAAT", "TTCAACAT", "TTCAAGAT", "TTCAATAA", "TTCAATAC",
        "TTCAATAG", "TTCAATAT", "TTCAATCT", "TTCAATGT", "TTCAATTT", "TTCACTAT", "TTCAGTAT", "TTCATTAT",
        "TTCCATAT", "TTCGATAT", "TTCTATAT", "TTGAATAT", "TTTAATAT"
    });
    EXPECT_EQ(
        NeighborsRecursive("TCCAATAT", 2),
        test2
    );

    EXPECT_EQ(
        NeighborsIterative("ACG", 1),
        test1
    );

    EXPECT_EQ(
        NeighborsIterative("TCCAATAT", 2),
        test2
    );
}

} // namespace