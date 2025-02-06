#pragma once

#include <algorithm>
#include <cstddef>
#include <ios>
#include <iosfwd>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <string>

namespace bs {

// keeps track of the first and last set bits
template <std::size_t N> class fl_bitset {
public:
  class reference;

  constexpr fl_bitset() noexcept {}

  constexpr fl_bitset(unsigned long long val) noexcept {
    std::size_t M = std::min(N, block_t_bitsize);
    data[0] |= val;
  }

  template <class charT = char, class traits = std::char_traits<charT>,
            class Allocator = std::allocator<charT>>

  explicit fl_bitset(
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

  // 20.9.2.2, bitset operations
  fl_bitset<N>& operator&=(const fl_bitset<N>& rhs) noexcept {
    for (std::size_t i = 0; i < num_blocks; i++) {
      data[i] &= rhs.data[i];
    }
    return *this;
  };

  fl_bitset<N>& operator|=(const fl_bitset<N>& rhs) noexcept {
    for (std::size_t i = 0; i < num_blocks; i++) {
      data[i] |= rhs.data[i];
    }
    return *this;
  };

  fl_bitset<N>& operator^=(const fl_bitset<N>& rhs) noexcept {
    for (std::size_t i = 0; i < num_blocks; i++) {
      data[i] ^= rhs.data[i];
    }
    return *this;
  };

  fl_bitset<N>& operator<<=(std::size_t pos) noexcept {
    namespace sv = std::ranges::views;
    for (std::size_t i : sv::iota(pos, N) | sv::reverse) {
      set_unchecked(i, (*this)[i - pos]);
    }
    for (std::size_t i = 0; i < pos; i++) {
      if (i % block_t_bitsize == 0 && i + block_t_bitsize < pos) {
        data[i / block_t_bitsize] = 0;
        i += block_t_bitsize - 1;
      } else {
        reset_unchecked(i);
      }
    }
    return *this;
  }

  fl_bitset<N>& operator>>=(std::size_t pos) noexcept {
    for (std::size_t i = 0; i + pos < N; i++) {
      set_unchecked(i, (*this)[i + pos]);
    }
    for (std::size_t i = N - pos; i < N; i++) {
      if (i % block_t_bitsize == 0) {
        data[i / block_t_bitsize] = 0;
        i += block_t_bitsize - 1;
      } else {
        reset_unchecked(i);
      }
    }
    return *this;
  };

  fl_bitset<N>& set() noexcept { return set_unchecked(); };

  fl_bitset<N>& set(std::size_t pos, bool val = true) {
    if (pos >= N)
      throw std::out_of_range{"Attempted to set bit out of range"};
    return set_unchecked(pos, val);
  }

  fl_bitset<N>& reset() noexcept { return reset_unchecked(); }

  fl_bitset<N>& reset(std::size_t pos) { return reset_unchecked(pos); }

  fl_bitset<N> operator~() const noexcept { return this->flip(); }

  fl_bitset<N>& flip() noexcept {
    for (std::size_t i = 0; i < num_blocks; i++) {
      data[i] = ~data[i];
    }
    return *this;
  }

  fl_bitset<N>& flip(std::size_t pos) {
    std::size_t block_idx = pos / block_t_bitsize;
    std::size_t bit = pos - (block_t_bitsize * block_idx);
    data[block_idx] ^= static_cast<block_t>(1) << bit;
    return *this;
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
    std::size_t sum = 0;
    for (std::size_t i = 0; i < num_blocks; i++) {
      sum += std::popcount(data[i]);
    }
    return sum;
  }

  constexpr std::size_t size() const noexcept { return N; }

  bool operator==(const fl_bitset<N>& rhs) const noexcept {
    for (std::size_t i = 0; i < num_blocks; i++) {
      if (data[i] != rhs.data[i])
        return false;
    }
    return true;
  }

  bool test(std::size_t pos) const {
    if (pos >= N)
      throw std::out_of_range{"Attempted to test bit out of range"};
    return (*this)[pos];
  }

  bool all() const noexcept {
    block_t mask = std::numeric_limits<block_t>::max();
    for (std::size_t i = 0; i < num_blocks; i++) {
      if (mask != (mask & data[i]))
        return false;
    }
    return true;
  }

  bool any() const noexcept {
    block_t mask = 0;
    for (std::size_t i = 0; i < num_blocks; i++) {
      mask |= data[i];
      if (mask)
        return true;
    }
    return false;
  }

  bool none() const noexcept { return !any(); }

  fl_bitset<N> operator<<(std::size_t pos) const noexcept {
    auto ret = *this;
    ret <<= pos;
    return ret;
  }

  fl_bitset<N> operator>>(std::size_t pos) const noexcept {
    auto ret = *this;
    ret >>= pos;
    return ret;
  }

  template <class charT, class traits>
  friend std::basic_istream<charT, traits>&
  operator>>(std::basic_istream<charT, traits>& is, fl_bitset<N>& x);

  // private:
  using block_t = std::size_t;
  constexpr static std::size_t block_t_bitsize = 8 * sizeof(block_t);
  constexpr static std::size_t num_blocks =
      (N + block_t_bitsize - 1) / block_t_bitsize;

  block_t data[num_blocks]{};

  fl_bitset<N>& set_unchecked();

  fl_bitset<N>& set_unchecked(std::size_t pos, bool val = true);

  fl_bitset<N>& reset_unchecked() noexcept;

  fl_bitset<N>& reset_unchecked(std::size_t pos);

  // need a non-noexcept impl to reuse in bitshift operators
};

template <std::size_t N> class fl_bitset<N>::reference {
public:
  friend fl_bitset;

  reference() = delete;
  reference(const reference&) = default;
  ~reference() = default;

  reference& operator=(bool x) noexcept {
    auto bit = static_cast<block_t>(1) << bit_idx;
    if (x) {
      block |= bit;
    } else {
      block &= std::numeric_limits<block_t>::max() ^ bit;
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

template <std::size_t N> fl_bitset<N>& fl_bitset<N>::set_unchecked() {
  block_t mask = std::numeric_limits<block_t>::max();
  for (std::size_t i = 0; i < num_blocks; i++) {
    data[i] |= mask;
  }
  return *this;
}

template <std::size_t N>
fl_bitset<N>& fl_bitset<N>::set_unchecked(std::size_t pos, bool val) {
  if (!val)
    return reset_unchecked(pos);
  std::size_t block_idx = pos / block_t_bitsize;
  std::size_t bit = pos - (block_t_bitsize * block_idx);
  data[block_idx] |= (static_cast<block_t>(1) << bit);
  return *this;
}

template <std::size_t N> fl_bitset<N>& fl_bitset<N>::reset_unchecked() noexcept {
  for (std::size_t i = 0; i < num_blocks; i++) {
    data[i] &= 0;
  }
  return *this;
}

template <std::size_t N>
fl_bitset<N>& fl_bitset<N>::reset_unchecked(std::size_t pos) {
  std::size_t block_idx = pos / block_t_bitsize;
  std::size_t bit = pos - (block_t_bitsize * block_idx);
  block_t mask = std::numeric_limits<block_t>::max();
  data[block_idx] &= mask ^ (static_cast<block_t>(1) << bit);
  return *this;
}

template <std::size_t N>
fl_bitset<N> operator&(const fl_bitset<N>& lhs, const fl_bitset<N>& rhs) noexcept {
  fl_bitset<N> res(lhs);
  res &= rhs;
  return res;
}

template <std::size_t N>
fl_bitset<N> operator|(const fl_bitset<N>& lhs, const fl_bitset<N>& rhs) noexcept {
  fl_bitset<N> res(lhs);
  res |= rhs;
  return res;
}

template <std::size_t N>
fl_bitset<N> operator^(const fl_bitset<N>& lhs, const fl_bitset<N>& rhs) noexcept {
  fl_bitset<N> res(lhs);
  res ^= rhs;
  return res;
};

template <class charT, class traits, std::size_t N>
std::basic_istream<charT, traits>&
operator>>(std::basic_istream<charT, traits>& is, fl_bitset<N>& x) {
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

  x = fl_bitset<N>{str};
  return is;
}

template <class charT, class traits, std::size_t N>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const fl_bitset<N>& x) {
  os << x.to_string();
  return os;
};

} // namespace bs

