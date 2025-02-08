// benchmark_count.cpp
#include "../include/bitset.hpp"
#include <benchmark/benchmark.h>
#include <execution>

constexpr std::size_t N = 1000000;

template<class T>
void do_it(benchmark::State& state, T& bits) {
  for (std::size_t i = 0; i < N / 2; ++i) {
    bits.set(i, true);
  }
  for (auto _ : state) {
    auto cnt = bits.count();
    benchmark::DoNotOptimize(cnt);
  }
}

// Variant 1: default
static void BM_Count_Default(benchmark::State& state) { 
  bs::bitset<N> bits; 
  do_it(state, bits);
}
BENCHMARK(BM_Count_Default);

// Variant 2: sequenced_policy
static void BM_Count_Seq(benchmark::State& state) {
  bs::bitset<N, std::execution::sequenced_policy, std::size_t> bits;
  do_it(state, bits);
}
BENCHMARK(BM_Count_Seq);

// Variant 3: parallel_policy
static void BM_Count_Par(benchmark::State& state) {
  bs::bitset<N, std::execution::parallel_policy, std::size_t> bits;
  do_it(state, bits);
}
BENCHMARK(BM_Count_Par);

// Variant 4: parallel_unsequenced_policy
static void BM_Count_ParUnseq(benchmark::State& state) {
  bs::bitset<N, std::execution::parallel_unsequenced_policy, std::size_t> bits;
  do_it(state, bits);
}
BENCHMARK(BM_Count_ParUnseq);

BENCHMARK_MAIN();
