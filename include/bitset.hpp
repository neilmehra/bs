#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <execution>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <ranges>
#include <stdexcept>
#include <string>

namespace bs {

template <class block_t>
concept C_BitOps = requires(block_t a, block_t b) {
  { block_t{0} } -> std::same_as<block_t>;
  { static_cast<block_t>(1) } -> std::same_as<block_t>;
  { a &= b } -> std::same_as<block_t&>;
  { a |= b } -> std::same_as<block_t&>;
  { a ^= b } -> std::same_as<block_t&>;
  { a | b } -> std::same_as<block_t>;
  { a& b } -> std::same_as<block_t>;
  { a ^ b } -> std::same_as<block_t>;
  { ~a } -> std::same_as<block_t>;
};

template <std::size_t N, class ExecutionPolicy, class block_t>
concept C_Bitset =
    (N > 0) && std::is_execution_policy_v<ExecutionPolicy> && C_BitOps<block_t>;

template <std::size_t N,
          class ExecutionPolicy = std::execution::sequenced_policy,
          class block_t = std::size_t>
  requires C_Bitset<N, ExecutionPolicy, block_t>
class bitset {
public:
  using backing_t = std::array<block_t, N>;
  using backing_it_t = backing_t::iterator;

  class reference;

  constexpr bitset() noexcept {}

  constexpr bitset(unsigned long long val) noexcept {
    std::size_t M = std::min(N, block_t_bitsize);
    data[0] |= val;
  }

  template <class charT = char, class traits = std::char_traits<charT>,
            class Allocator = std::allocator<charT>>

  explicit bitset(
      const std::basic_string<charT, traits, Allocator>& str,
      typename std::basic_string<charT, traits, Allocator>::size_type pos = 0,
      typename std::basic_string<charT, traits, Allocator>::size_type n =
          std::basic_string<charT, traits, Allocator>::npos,
      charT zero = charT('0'), charT one = charT('1')) {

    using str_t = std::basic_string<charT, traits, Allocator>;
    if (pos > str.size())
      throw std::out_of_range{"string starting index out of range"};

    std::size_t idx = 0;
    while (idx < std::min(N, n == str_t::npos ? str.size() : n)) {
      if (str[pos + idx] != zero && str[pos + idx] != one) {
        throw std::invalid_argument{"non-0/1 char in str"};
      } else if (str[pos + idx] == one) {
        set(pos + idx);
      }
      idx++;
    }
  }

  backing_it_t begin() noexcept { return data.begin(); }
  backing_it_t end() noexcept { return data.end(); }

  backing_it_t begin() const noexcept { return data.begin(); }
  backing_it_t end() const noexcept { return data.end(); }

  template <class UnaryOp> bitset& apply_op(UnaryOp f) {
    std::for_each(execution_policy, data.begin(), data.end(), f);
    return *this;
  }

  template <class UnaryOp>
  bitset& apply_op(backing_it_t begin, backing_it_t end, UnaryOp f) {
    std::for_each(execution_policy, begin, end, f);
    return *this;
  }

  template <class BinaryOp> bitset& apply_zip(bitset& rhs, BinaryOp f) {
    auto z = std::views::zip(data, rhs.data);
    std::for_each(execution_policy, z.begin(), z.end(), f);
    return *this;
  }

  template <class BinaryOp>
  bitset& apply_zip(bitset& rhs, backing_it_t begin, backing_it_t end,
                    BinaryOp f) {
    auto beg_dist = std::distance(data.begin(), begin);
    auto end_dist = std::distance(data.begin(), end);
    auto z = std::views::zip(data, rhs.data);
    std::for_each(execution_policy, z.begin() + beg_dist, z.begin() + end_dist,
                  f);
    return *this;
  }

  // 20.9.2.2, bitset operations
  bitset& operator&=(const bitset& rhs) noexcept {
    return apply_zip(const_cast<bitset&>(rhs), [](const auto& v) {
      auto& [lhs, rhs] = v;
      lhs &= rhs;
    });
  };

  bitset& operator|=(const bitset& rhs) noexcept {
    return apply_zip(const_cast<bitset&>(rhs), [](const auto& v) {
      auto& [lhs, rhs] = v;
      lhs |= rhs;
    });
  };

  bitset& operator^=(const bitset& rhs) noexcept {
    return apply_zip(const_cast<bitset&>(rhs), [](const auto& v) {
      auto& [lhs, rhs] = v;
      lhs ^= rhs;
    });
  };

  bitset& operator<<=(std::size_t pos) noexcept {
    namespace sv = std::ranges::views;
    for (std::size_t i : sv::iota(pos, N) | sv::reverse) {
      set_unchecked(i, (*this)[i - pos]);
    }

    std::size_t pos_block_idx = pos / block_t_bitsize;
    reset(data.begin(), data.begin() + pos_block_idx);

    std::size_t bit_idx = pos - (pos_block_idx * block_t_bitsize);
    block_t mask = ~((static_cast<block_t>(1) << bit_idx) - 1);
    data[pos_block_idx] &= mask;

    return *this;
  }

  bitset& operator>>=(std::size_t pos) noexcept {
    for (std::size_t i = 0; i + pos < N; i++) {
      set_unchecked(i, (*this)[i + pos]);
    }

    std::size_t pos_block_idx = pos / block_t_bitsize;
    reset(data.begin() + pos_block_idx + 1, data.end());

    std::size_t bit_idx = pos - (pos_block_idx * block_t_bitsize);

    block_t mask =
        (static_cast<block_t>(1) << (std::min(block_t_bitsize, N) - bit_idx)) -
        1;
    data[pos_block_idx] &= mask;

    return *this;
  };

  bitset& set() noexcept { return set(data.begin(), data.end()); };

  bitset& set(std::size_t pos, bool val = true) {
    if (pos >= N)
      throw std::out_of_range{"Attempted to set bit out of range"};
    return set_unchecked(pos, val);
  }

  bitset& set(backing_it_t begin, backing_it_t end) {
    block_t mask = ~block_t{0};
    return this->apply_op(begin, end, [mask](auto& v) { v = mask; });
  }

  bitset& reset() noexcept { return reset(data.begin(), data.end()); }

  bitset& reset(std::size_t pos) {
    if (pos >= N) {
      throw std::out_of_range{"Attempted to reset bit out of range"};
    }
    return reset_unchecked(pos);
  }

  bitset& reset(backing_it_t begin, backing_it_t end) {
    return apply_op(begin, end, [](auto& v) { v = 0; });
  }

  bitset operator~() const noexcept { return this->flip(); }

  bitset& flip() noexcept { return flip(data.begin(), data.end()); }

  bitset& flip(std::size_t pos) {
    std::size_t block_idx = pos / block_t_bitsize;
    std::size_t bit = pos - (block_t_bitsize * block_idx);
    data[block_idx] ^= static_cast<block_t>(1) << bit;
    return *this;
  }

  bitset& flip(backing_it_t first, backing_it_t last) {
    return apply_op(first, last, [](auto& v) { v = ~v; });
  }

  // element access
  constexpr bool operator[](std::size_t pos) const {
    std::size_t block_idx = pos / block_t_bitsize;
    std::size_t bit = pos - (block_t_bitsize * block_idx);
    return (data[block_idx] & (static_cast<block_t>(1) << bit)) != 0;
  }

  // for b[i];
  reference operator[](std::size_t pos) {
    std::size_t block_idx = pos / block_t_bitsize;
    std::size_t bit = pos - (block_t_bitsize * block_idx);
    return reference{data[block_idx], bit};
  }

  // for b[i];
  unsigned long to_ulong() const {
    for (std::size_t i = 1; i < num_blocks; i++) {
      if (data[i])
        throw std::overflow_error{"Incurred overflow upon attempting to "
                                  "convert bitset to unsigned long"};
    }
    return static_cast<unsigned long>(data[0]);
  }

  unsigned long long to_ullong() const { return to_ulong(); }

  template <class charT = char, class traits = std::char_traits<charT>,
            class Allocator = std::allocator<charT>>
  std::basic_string<charT, traits, Allocator>
  to_string(charT zero = charT('0'), charT one = charT('1')) const {
    std::basic_string<charT, traits, Allocator> ret(N, zero);

    for (std::size_t i = 0; i < N; i++) {
      if (test(i))
        ret[N - i - 1] = one;
    }
    return ret;
  }

  std::size_t count() const noexcept {
    return std::reduce(execution_policy, data.begin(), data.end(),
                       std::size_t{0}, [](const auto& accum, const auto& v) {
                         return accum + std::popcount(v);
                       });
  }

  constexpr std::size_t size() const noexcept { return N; }

  bool operator==(const bitset& rhs) const noexcept {
    auto z = std::views::zip(this->data, rhs.data);
    return std::all_of(execution_policy, z.begin(), z.end(), [](const auto& v) {
      const auto& [lhs, rhs] = v;
      return lhs == rhs;
    });
  }

  bool test(std::size_t pos) const {
    if (pos >= N)
      throw std::out_of_range{"Attempted to test bit out of range"};
    return (*this)[pos];
  }

  bool all() const noexcept {
    block_t mask = ~(block_t{0});
    return std::all_of(execution_policy, data.begin(), data.end(),
                       [mask](const auto& v) { return v == mask; });
  }

  bool any() const noexcept {
    return std::any_of(execution_policy, data.begin(), data.end(),
                       [](const auto& v) { return v != 0; });
  }

  bool none() const noexcept { return !any(); }

  bitset operator<<(std::size_t pos) const noexcept {
    auto ret = *this;
    ret <<= pos;
    return ret;
  }

  bitset operator>>(std::size_t pos) const noexcept {
    auto ret = *this;
    ret >>= pos;
    return ret;
  }

  template <class charT, class traits>
  friend std::basic_istream<charT, traits>&
  operator>>(std::basic_istream<charT, traits>& is, bitset& x);

private:
  constexpr static ExecutionPolicy execution_policy{};
  constexpr static std::size_t block_t_bitsize = 8 * sizeof(block_t);
  constexpr static std::size_t num_blocks =
      (N + block_t_bitsize - 1) / block_t_bitsize;

  backing_t data{};

  // need a non-noexcept impl to reuse in bitshift operators
  bitset& set_unchecked(std::size_t pos, bool val) {
    if (!val)
      return reset_unchecked(pos);
    std::size_t block_idx = pos / block_t_bitsize;
    std::size_t bit = pos - (block_t_bitsize * block_idx);
    data[block_idx] |= (static_cast<block_t>(1) << bit);
    return *this;
  }

  bitset& reset_unchecked(std::size_t pos) {
    std::size_t block_idx = pos / block_t_bitsize;
    std::size_t bit = pos - (block_t_bitsize * block_idx);
    block_t mask = ~block_t{0};
    data[block_idx] &= mask ^ (static_cast<block_t>(1) << bit);
    return *this;
  }
};

template <std::size_t N, class ExecutionPolicy, class block_t>
  requires C_Bitset<N, ExecutionPolicy, block_t>
class bitset<N, ExecutionPolicy, block_t>::reference {
public:
  friend class bitset;

  reference() = delete;
  reference(const reference&) = default;
  ~reference() = default;

  reference& operator=(bool x) noexcept {
    auto bit = static_cast<block_t>(1) << bit_idx;
    if (x) {
      block |= bit;
    } else {
      block &= (~block_t{0}) ^ bit;
    }
    return *this;
  };

  reference& operator=(const reference& rhs) noexcept {
    block = rhs.block;
    bit_idx = rhs.bit_idx;
    return *this;
  };

  bool operator~() const noexcept { return !static_cast<bool>(*this); }

  operator bool() const noexcept {
    return (block & (static_cast<block_t>(1) << bit_idx)) != 0;
  };

  reference& flip() noexcept {
    *this = ~(*this);
    return *this;
  };

private:
  block_t& block;
  std::size_t bit_idx;

  reference(block_t& block_, std::size_t bit_idx_) noexcept
      : block(block_), bit_idx(bit_idx_) {}
};

template <std::size_t N, class ExecutionPolicy, class block_t>
bitset<N, ExecutionPolicy, block_t>
operator&(const bitset<N, ExecutionPolicy, block_t>& lhs,
          const bitset<N, ExecutionPolicy, block_t>& rhs) noexcept {
  bitset<N, ExecutionPolicy, block_t> res(lhs);
  res &= rhs;
  return res;
}

template <std::size_t N, class ExecutionPolicy, class block_t>
bitset<N, ExecutionPolicy, block_t>
operator|(const bitset<N, ExecutionPolicy, block_t>& lhs,
          const bitset<N, ExecutionPolicy, block_t>& rhs) noexcept {
  bitset<N, ExecutionPolicy, block_t> res(lhs);
  res |= rhs;
  return res;
}

template <std::size_t N, class ExecutionPolicy, class block_t>
bitset<N, ExecutionPolicy, block_t>
operator^(const bitset<N, ExecutionPolicy, block_t>& lhs,
          const bitset<N, ExecutionPolicy, block_t>& rhs) noexcept {
  bitset<N, ExecutionPolicy, block_t> res(lhs);
  res ^= rhs;
  return res;
};

template <class charT, class traits, std::size_t N, class ExecutionPolicy,
          class block_t>
std::basic_istream<charT, traits>&
operator>>(std::basic_istream<charT, traits>& is,
           bitset<N, ExecutionPolicy, block_t>& x) {
  // basic_istream::widen => Converts a character to its equivalent in the
  // current locale. The result is converted from char to character type used
  // within the stream if needed.
  auto zero = is.widen('0'), one = is.widen('1'), ws = is.widen(' ');

  std::size_t idx = 0;
  std::basic_string<charT, traits> str(N, zero);
  charT test_c;

  while (is >> test_c && test_c == ws) {
  }; // clang-format moment

  while (str.size() < N && is >> test_c) {
    if (test_c == zero || test_c == one) {
      str[idx++] = test_c;
    } else {
      break;
    }
  }

  if (N > 0 && idx == 0) {
    is.setstate(std::ios_base::failbit);
  }

  x = bitset<N, ExecutionPolicy, block_t>{str};
  return is;
}

template <class charT, class traits, std::size_t N, class ExecutionPolicy,
          class block_t>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os,
           const bitset<N, ExecutionPolicy, block_t>& x) {
  os << x.to_string();
  return os;
};

} // namespace bs
