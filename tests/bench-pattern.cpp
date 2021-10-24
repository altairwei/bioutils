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
    Benchmark for PatternCount
*/

template <typename ...ExtraArgs>
void BenchPatternCount(benchmark::State& state, ExtraArgs&&... extra_args) {
    auto sequence_length = state.range(0);
    auto pattern_length = 10;

    std::string genome = random_sequence(sequence_length);
    std::string pattern = random_sequence(pattern_length);

    for (auto _ : state) {
        PatternCount(genome, pattern, extra_args...);
    }
}

BENCHMARK_CAPTURE(BenchPatternCount, Fast, AlgorithmEfficiency::Fast)->RangeMultiplier(2)->Range(1024, 1024<<12);
BENCHMARK_CAPTURE(BenchPatternCount, Fastest, AlgorithmEfficiency::Fastest)->RangeMultiplier(2)->Range(1024, 1024<<12);

/*
    Benchmark for FindClumps
*/

template <typename ...ExtraArgs>
void BenchFindClumps(benchmark::State& state, ExtraArgs&&... extra_args) {
    auto sequence_length = state.range(0);
    std::string genome = random_sequence(sequence_length);

    for (auto _ : state) {
        FindClumps(genome, 10, 500, 20, extra_args...);
    }
}

//BENCHMARK_CAPTURE(BenchFindClumps, Slow, 10, 500, 20, AlgorithmEfficiency::Slow)->RangeMultiplier(2)->Range(1024, 1024<<4);
//BENCHMARK_CAPTURE(BenchFindClumps, Fast, AlgorithmEfficiency::Fast)->RangeMultiplier(2)->Range(1024, 1024<<10);
BENCHMARK_CAPTURE(BenchFindClumps, Faster, AlgorithmEfficiency::Faster)->RangeMultiplier(2)->Range(1024, 1024<<12);
BENCHMARK_CAPTURE(BenchFindClumps, Fastest, AlgorithmEfficiency::Fastest)->RangeMultiplier(2)->Range(1024, 1024<<12);