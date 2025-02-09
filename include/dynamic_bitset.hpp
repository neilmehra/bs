#pragma once

#include <algorithm>
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

template <class ExecutionPolicy, class block_t>
concept C_Bitset =
    std::is_execution_policy_v<ExecutionPolicy> && C_BitOps<block_t>;

template <class ExecutionPolicy = std::execution::sequenced_policy,
          class block_t = std::size_t>
  requires C_Bitset<ExecutionPolicy, block_t>
class dbitset {
public:
  using backing_t = std::vector<block_t>;
  using backing_it_t = backing_t::iterator;

  class reference;

  dbitset() noexcept {}

  dbitset(unsigned long long val) {
    data.reserve(1);
    data.push_back(val);
  }

  template <class charT = char, class traits = std::char_traits<charT>,
            class Allocator = std::allocator<charT>>
  explicit dbitset(
      const std::basic_string<charT, traits, Allocator>& str,
      typename std::basic_string<charT, traits, Allocator>::size_type pos = 0,
      typename std::basic_string<charT, traits, Allocator>::size_type n =
          std::basic_string<charT, traits, Allocator>::npos,
      charT zero = charT('0'), charT one = charT('1'));

  backing_it_t begin() noexcept { return data.begin(); }
  backing_it_t end() noexcept { return data.end(); }

  backing_it_t begin() const noexcept { return data.begin(); }
  backing_it_t end() const noexcept { return data.end(); }

  template <class UnaryOp> dbitset& apply_op(UnaryOp f);

  template <class UnaryOp>
  dbitset& apply_op(backing_it_t begin, backing_it_t end, UnaryOp f);

  template <class BinaryOp> dbitset& apply_zip(dbitset& rhs, BinaryOp f);

  template <class BinaryOp>
  dbitset& apply_zip(dbitset& rhs, backing_it_t begin, backing_it_t end,
                     BinaryOp f);

  dbitset& operator&=(const dbitset& rhs) noexcept {
    return apply_zip(const_cast<dbitset&>(rhs), [](const auto& v) {
      auto& [lhs, rhs] = v;
      lhs &= rhs;
    });
  }

  dbitset& operator|=(const dbitset& rhs) noexcept {
    return apply_zip(const_cast<dbitset&>(rhs), [](const auto& v) {
      auto& [lhs, rhs] = v;
      lhs |= rhs;
    });
  }

  dbitset& operator^=(const dbitset& rhs) noexcept {
    return apply_zip(const_cast<dbitset&>(rhs), [](const auto& v) {
      auto& [lhs, rhs] = v;
      lhs ^= rhs;
    });
  }

  dbitset& operator<<=(std::size_t pos) noexcept;

  dbitset& operator>>=(std::size_t pos) noexcept;

  dbitset operator~() const noexcept { return this->flip(); }

  constexpr bool operator[](std::size_t pos) const;

  reference operator[](std::size_t pos);

  bool operator==(const dbitset& rhs) const noexcept;

  dbitset operator<<(std::size_t pos) const noexcept {
    auto ret = *this;
    return ret <<= pos;
  }

  dbitset operator>>(std::size_t pos) const noexcept {
    auto ret = *this;
    return ret >>= pos;
  }

  dbitset& set() noexcept { return set(data.begin(), data.end()); }

  dbitset& set(std::size_t pos, bool val = true) {
    if (pos >= size())
      throw std::out_of_range{"Attempted to set bit out of range"};
    return set_unchecked(pos, val);
  }

  dbitset& set(backing_it_t begin, backing_it_t end) {
    block_t mask = ~block_t{0};
    return this->apply_op(begin, end, [mask](auto& v) { v = mask; });
  }

  dbitset& reset() noexcept { return reset(data.begin(), data.end()); }

  dbitset& reset(std::size_t pos) {

    if (pos >= size()) {
      throw std::out_of_range{"Attempted to reset bit out of range"};
    }
    return reset_unchecked(pos);
  }

  dbitset& reset(backing_it_t begin, backing_it_t end) {

    return apply_op(begin, end, [](auto& v) { v = 0; });
  }

  dbitset& flip() noexcept { return flip(data.begin(), data.end()); }

  dbitset& flip(std::size_t pos);

  dbitset& flip(backing_it_t first, backing_it_t last) {
    return apply_op(first, last, [](auto& v) { v = ~v; });
  }

  // data access

  std::size_t count() const noexcept;

  std::size_t size() const noexcept { return size_; }

  bool test(std::size_t pos) const {

    if (pos >= size())
      throw std::out_of_range{"Attempted to test bit out of range"};
    return (*this)[pos];
  }

  bool all() const noexcept;

  bool any() const noexcept {
    return std::any_of(execution_policy, data.begin(), data.end(),
                       [](const auto& v) { return v != 0; });
  }

  bool none() const noexcept { return !any(); }

  unsigned long to_ulong() const;

  unsigned long long to_ullong() const { return to_ulong(); }

  template <class charT = char, class traits = std::char_traits<charT>,
            class Allocator = std::allocator<charT>>
  std::basic_string<charT, traits, Allocator>
  to_string(charT zero = charT('0'), charT one = charT('1')) const;

  template <class charT, class traits>
  friend std::basic_istream<charT, traits>&
  operator>>(std::basic_istream<charT, traits>& is, dbitset& x);

  // dynamic methods here

  void push_back(bool b) {

    reserve();
    set(size_++, b);
  }

private:
  constexpr static ExecutionPolicy execution_policy{};
  constexpr static std::size_t block_t_bitsize = 8 * sizeof(block_t);

  std::size_t size_; // for num bits
  backing_t data{};

  void reserve() {
    if ((size_ + 1) % block_t_bitsize == 0) {
      data.push_back(0);
    }
  }

  // need a non-noexcept impl to reuse in bitshift operators
  dbitset& set_unchecked(std::size_t pos, bool val);

  dbitset& reset_unchecked(std::size_t pos);
};

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
template <class charT, class traits, class Allocator>
dbitset<ExecutionPolicy, block_t>::dbitset(
    const std::basic_string<charT, traits, Allocator>& str,
    typename std::basic_string<charT, traits, Allocator>::size_type pos,
    typename std::basic_string<charT, traits, Allocator>::size_type n,
    charT zero, charT one) {

  using str_t = std::basic_string<charT, traits, Allocator>;
  if (pos > str.size())
    throw std::out_of_range{"string starting index out of range"};

  std::size_t num_chars = n == str_t::npos ? str.size() : n;
  data.reserve(num_chars);

  std::size_t idx = 0;
  while (idx < num_chars) {
    if (str[pos + idx] != zero && str[pos + idx] != one) {
      throw std::invalid_argument{"non-0/1 char in str"};
    } else if (str[pos + idx] == one) {
      set(pos + idx);
    }
    idx++;
  }
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
template <class UnaryOp>
dbitset<ExecutionPolicy, block_t>&
dbitset<ExecutionPolicy, block_t>::apply_op(UnaryOp f) {
  std::for_each(execution_policy, data.begin(), data.end(), f);
  return *this;
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
template <class UnaryOp>
dbitset<ExecutionPolicy, block_t>&
dbitset<ExecutionPolicy, block_t>::apply_op(backing_it_t begin,
                                            backing_it_t end, UnaryOp f) {
  std::for_each(execution_policy, begin, end, f);
  return *this;
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
template <class BinaryOp>
dbitset<ExecutionPolicy, block_t>&
dbitset<ExecutionPolicy, block_t>::apply_zip(dbitset& rhs, BinaryOp f) {
  auto z = std::views::zip(data, rhs.data);
  std::for_each(execution_policy, z.begin(), z.end(), f);
  return *this;
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
template <class BinaryOp>
dbitset<ExecutionPolicy, block_t>&
dbitset<ExecutionPolicy, block_t>::apply_zip(dbitset& rhs, backing_it_t begin,
                                             backing_it_t end, BinaryOp f) {
  auto beg_dist = std::distance(data.begin(), begin);
  auto end_dist = std::distance(data.begin(), end);
  auto z = std::views::zip(data, rhs.data);
  std::for_each(execution_policy, z.begin() + beg_dist, z.begin() + end_dist,
                f);
  return *this;
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
dbitset<ExecutionPolicy, block_t>&
dbitset<ExecutionPolicy, block_t>::operator<<=(std::size_t pos) noexcept {
  namespace sv = std::ranges::views;
  for (std::size_t i : sv::iota(pos, size()) | sv::reverse) {
    set_unchecked(i, (*this)[i - pos]);
  }

  std::size_t pos_block_idx = pos / block_t_bitsize;
  reset(data.begin(), data.begin() + pos_block_idx);

  std::size_t bit_idx = pos - (pos_block_idx * block_t_bitsize);
  block_t mask = ~((static_cast<block_t>(1) << bit_idx) - 1);
  data[pos_block_idx] &= mask;

  return *this;
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
dbitset<ExecutionPolicy, block_t>&
dbitset<ExecutionPolicy, block_t>::operator>>=(std::size_t pos) noexcept {
  for (std::size_t i = 0; i + pos < size(); i++) {
    set_unchecked(i, (*this)[i + pos]);
  }

  std::size_t pos_block_idx = pos / block_t_bitsize;
  reset(data.begin() + pos_block_idx + 1, data.end());

  std::size_t bit_idx = pos - (pos_block_idx * block_t_bitsize);

  block_t mask = (static_cast<block_t>(1)
                  << (std::min(block_t_bitsize, size()) - bit_idx)) -
                 1;
  data[pos_block_idx] &= mask;

  return *this;
};

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
dbitset<ExecutionPolicy, block_t>&
dbitset<ExecutionPolicy, block_t>::flip(std::size_t pos) {
  std::size_t block_idx = pos / block_t_bitsize;
  std::size_t bit = pos - (block_t_bitsize * block_idx);
  data[block_idx] ^= static_cast<block_t>(1) << bit;
  return *this;
}

// element access
template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
constexpr bool
dbitset<ExecutionPolicy, block_t>::operator[](std::size_t pos) const {
  std::size_t block_idx = pos / block_t_bitsize;
  std::size_t bit = pos - (block_t_bitsize * block_idx);
  return (data[block_idx] & (static_cast<block_t>(1) << bit)) != 0;
}

// for b[i];
template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
dbitset<ExecutionPolicy, block_t>::reference
dbitset<ExecutionPolicy, block_t>::operator[](std::size_t pos) {
  std::size_t block_idx = pos / block_t_bitsize;
  std::size_t bit = pos - (block_t_bitsize * block_idx);
  return reference{data[block_idx], bit};
}

// for b[i];
template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
unsigned long dbitset<ExecutionPolicy, block_t>::to_ulong() const {
  for (std::size_t i = 1; i < data.size(); i++) {
    if (data[i])
      throw std::overflow_error{"Incurred overflow upon attempting to "
                                "convert dbitset to unsigned long"};
  }
  return static_cast<unsigned long>(data[0]);
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
template <class charT, class traits, class Allocator>
std::basic_string<charT, traits, Allocator>
dbitset<ExecutionPolicy, block_t>::to_string(charT zero, charT one) const {
  std::basic_string<charT, traits, Allocator> ret(size(), zero);

  for (std::size_t i = 0; i < size(); i++) {
    if (test(i))
      ret[size() - i - 1] = one;
  }
  return ret;
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
std::size_t dbitset<ExecutionPolicy, block_t>::count() const noexcept {
  return std::reduce(execution_policy, data.begin(), data.end(), std::size_t{0},
                     [](const auto& accum, const auto& v) {
                       return accum + std::popcount(v);
                     });
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
bool dbitset<ExecutionPolicy, block_t>::operator==(
    const dbitset& rhs) const noexcept {
  auto z = std::views::zip(this->data, rhs.data);
  return std::all_of(execution_policy, z.begin(), z.end(), [](const auto& v) {
    const auto& [lhs, rhs] = v;
    return lhs == rhs;
  });
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
bool dbitset<ExecutionPolicy, block_t>::all() const noexcept {
  block_t mask = ~(block_t{0});
  return std::all_of(execution_policy, data.begin(), data.end(),
                     [mask](const auto& v) { return v == mask; });
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
dbitset<ExecutionPolicy, block_t>&
dbitset<ExecutionPolicy, block_t>::set_unchecked(std::size_t pos, bool val) {
  if (!val)
    return reset_unchecked(pos);
  std::size_t block_idx = pos / block_t_bitsize;
  std::size_t bit = pos - (block_t_bitsize * block_idx);
  data[block_idx] |= (static_cast<block_t>(1) << bit);
  return *this;
}

template <class ExecutionPolicy, class block_t>
  requires C_Bitset<ExecutionPolicy, block_t>
dbitset<ExecutionPolicy, block_t>&
dbitset<ExecutionPolicy, block_t>::reset_unchecked(std::size_t pos) {
  std::size_t block_idx = pos / block_t_bitsize;
  std::size_t bit = pos - (block_t_bitsize * block_idx);
  block_t mask = ~block_t{0};
  data[block_idx] &= mask ^ (static_cast<block_t>(1) << bit);
  return *this;
}

template <class ExecutionPolicy, class block_t>
dbitset<ExecutionPolicy, block_t>
operator&(const dbitset<ExecutionPolicy, block_t>& lhs,
          const dbitset<ExecutionPolicy, block_t>& rhs) noexcept {
  dbitset<ExecutionPolicy, block_t> res(lhs);
  res &= rhs;
  return res;
}

template <class ExecutionPolicy, class block_t>
dbitset<ExecutionPolicy, block_t>
operator|(const dbitset<ExecutionPolicy, block_t>& lhs,
          const dbitset<ExecutionPolicy, block_t>& rhs) noexcept {
  dbitset<ExecutionPolicy, block_t> res(lhs);
  res |= rhs;
  return res;
}

template <class ExecutionPolicy, class block_t>
dbitset<ExecutionPolicy, block_t>
operator^(const dbitset<ExecutionPolicy, block_t>& lhs,
          const dbitset<ExecutionPolicy, block_t>& rhs) noexcept {
  dbitset<ExecutionPolicy, block_t> res(lhs);
  res ^= rhs;
  return res;
};

template <class charT, class traits, std::size_t N, class ExecutionPolicy,
          class block_t>
std::basic_istream<charT, traits>&
operator>>(std::basic_istream<charT, traits>& is,
           dbitset<ExecutionPolicy, block_t>& x) {
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

  x = dbitset<ExecutionPolicy, block_t>{str};
  return is;
}

template <class charT, class traits, class ExecutionPolicy, class block_t>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os,
           const dbitset<ExecutionPolicy, block_t>& x) {
  os << x.to_string();
  return os;
};

} // namespace bs
