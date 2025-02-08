// benchmark_apply_op.cpp
#include "../include/bitset.hpp"
#include <benchmark/benchmark.h>
#include <execution>

constexpr std::size_t N = 1000000;

template <class T> void do_it(benchmark::State& state, T& bits) {
  for (std::size_t i = 0; i < N / 2; ++i) {
    bits.set(i, true);
  }
  for (auto _ : state) {
    auto& res = bits.apply_op([](auto& v) { ++v; });
    benchmark::DoNotOptimize(res);
  }
}

// Variant 1: default
static void BM_ApplyOp_Default(benchmark::State& state) {
  bs::bitset<N> bits;
  do_it(state, bits);
}
BENCHMARK(BM_ApplyOp_Default);

// Variant 2: sequenced_policy
static void BM_ApplyOp_Seq(benchmark::State& state) {
  bs::bitset<N, std::execution::sequenced_policy, std::size_t> bits;
  do_it(state, bits);
}
BENCHMARK(BM_ApplyOp_Seq);

// Variant 3: parallel_policy
static void BM_ApplyOp_Par(benchmark::State& state) {
  bs::bitset<N, std::execution::parallel_policy, std::size_t> bits;
  do_it(state, bits);
}
BENCHMARK(BM_ApplyOp_Par);

// Variant 4: parallel_unsequenced_policy
static void BM_ApplyOp_ParUnseq(benchmark::State& state) {
  bs::bitset<N, std::execution::parallel_unsequenced_policy, std::size_t> bits;
  do_it(state, bits);
}
BENCHMARK(BM_ApplyOp_ParUnseq);

BENCHMARK_MAIN();
