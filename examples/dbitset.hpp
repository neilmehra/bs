#pragma once

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

namespace dbitset {

template <class Block, class Allocator> class dbitset {
public:
  using block_t = Block;
  using backing_t = std::vector<Block, Allocator>;

  dbitset();
  dbitset(unsigned long);
  dbitset(unsigned long long) noexcept;

  [[nodiscard]] bool operator==(const dbitset&) const;
  [[nodiscard]] bool operator!=(const dbitset&) const;

  // element access
  [[nodiscard]] bool operator[](std::size_t) const;
  [[nodiscard]] bool test(std::size_t) const;
  [[nodiscard]] bool any() const noexcept;
  [[nodiscard]] bool none() const noexcept;
  [[nodiscard]] bool all() const noexcept;
  [[nodiscard]] std::size_t count() const noexcept;

  [[nodiscard]] bool front() const noexcept;
  [[nodiscard]] bool back() const noexcept;
  [[nodiscard]] Block& front_block() const noexcept;
  [[nodiscard]] Block& back_block() const noexcept;

  // capacity
  [[nodiscard]] std::size_t size() const noexcept;

  dbitset& operator&=(const dbitset&) noexcept;
  dbitset& operator|=(const dbitset&) noexcept;
  dbitset& operator^=(const dbitset&) noexcept;
  dbitset& operator~() noexcept;

  dbitset& operator<<=(std::size_t) noexcept;
  dbitset& operator>>=(std::size_t) noexcept;
  dbitset& operator<<(std::size_t) noexcept;
  dbitset& operator>>(std::size_t) noexcept;

  dbitset& set() noexcept;
  dbitset& set(std::size_t);

  dbitset& reset() noexcept;
  dbitset& reset(std::size_t);

  dbitset& flip() noexcept;
  dbitset& flip(std::size_t);

  // conversions
  [[nodiscard]] std::string to_string() const noexcept;
  [[nodiscard]] unsigned long to_ulong() const;
  [[nodiscard]] unsigned long long to_ullong() const;

private:
  std::size_t size_;
};

} // namespace dbitset
