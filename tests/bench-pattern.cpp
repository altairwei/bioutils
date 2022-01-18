#include <random>
#include <iterator>
#include <algorithm>

#include "benchmark/benchmark.h"

#include "pattern.h"

using namespace bioutils::algorithms;

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

static const std::array<char, 4> NUCLEOTIDES = {'A', 'C', 'G', 'T'};

template <typename T>
static std::string random_sequence(const T length) {
    std::string seq(length, 'A');
    std::generate(seq.begin(), seq.end(), [] {
        return *select_randomly(NUCLEOTIDES.begin(), NUCLEOTIDES.end());
    });
    return seq;
}

/*
 * Benchmark for PatternCount
 * ——————————————————————————————————————————————————
 */

typedef std::size_t (*PatternCountFuncPtr)(std::string_view, std::string_view);
void BenchPatternCount(benchmark::State& state, PatternCountFuncPtr fun) {
    auto sequence_length = state.range(0);
    auto pattern_length = 10;

    std::string genome = random_sequence(sequence_length);
    std::string pattern = random_sequence(pattern_length);

    for (auto _ : state) {
        fun(genome, pattern);
    }
}

BENCHMARK_CAPTURE(BenchPatternCount, BruteForce, PatternCount_BF)->RangeMultiplier(2)->Range(1024, 1024<<12);
BENCHMARK_CAPTURE(BenchPatternCount, RabinKarp, PatternCount_RK)->RangeMultiplier(2)->Range(1024, 1024<<12);

/*
 * Benchmark for FrequentWords
 * ——————————————————————————————————————————————————
 */

typedef std::set<std::string> (*FrequentWordsFuncPtr)(std::string_view, int);
void BenchFrequentWords(benchmark::State& state, FrequentWordsFuncPtr fun) {
    auto sequence_length = state.range(0);
    std::string genome = random_sequence(sequence_length);

    for (auto _ : state) {
        fun(genome, 6);
    }
}

BENCHMARK_CAPTURE(BenchFrequentWords, ByPerfectHash, FrequentWordsByPerfectHash)->RangeMultiplier(2)->Range(1024, 1024<<12);
BENCHMARK_CAPTURE(BenchFrequentWords, BySorting, FrequentWordsBySorting)->RangeMultiplier(2)->Range(1024, 1024<<12);
BENCHMARK_CAPTURE(BenchFrequentWords, ByStdHash, FrequentWordsByStdHash)->RangeMultiplier(2)->Range(1024, 1024<<12);

/*
 * Benchmark for FindClumps
 * ——————————————————————————————————————————————————
 */

typedef std::set<std::string> (*FindClumpsFuncPtr)(std::string_view, int, int, int);
void BenchFindClumps(benchmark::State& state, FindClumpsFuncPtr fun) {
    auto sequence_length = state.range(0);
    std::string genome = random_sequence(sequence_length);

    for (auto _ : state) {
        fun(genome, 10, 500, 20);
    }
}

BENCHMARK_CAPTURE(BenchFindClumps, WithPerfectHash, FindClumpsBetterWithPerfectHash)->RangeMultiplier(2)->Range(1024, 1024<<12);
BENCHMARK_CAPTURE(BenchFindClumps, WithStdHash, FindClumpsBetterWithStdHash)->RangeMultiplier(2)->Range(1024, 1024<<12);