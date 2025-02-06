#include "../include/bitset.hpp"
#include <gtest/gtest.h>
#include <string>

TEST(BitsetTest, DefaultConstructor) {
  bs::bitset<8> b;
  for (int i = 0; i < 8; i++) {
    EXPECT_FALSE(b[i]);
  }
}

TEST(BitsetTest, ULongLongConstructor) {
  bs::bitset<8> b(0b10101010ULL);
  EXPECT_TRUE(b[1]);
  EXPECT_TRUE(b[3]);
  EXPECT_TRUE(b[5]);
  EXPECT_TRUE(b[7]);
  EXPECT_FALSE(b[0]);
  EXPECT_FALSE(b[2]);
  EXPECT_FALSE(b[4]);
  EXPECT_FALSE(b[6]);
}

TEST(BitsetTest, StringConstructor) {
  std::string s = "10110010";
  bs::bitset<8> b(s);
  EXPECT_TRUE(b[0]);
  EXPECT_TRUE(b[2]);
  EXPECT_TRUE(b[3]);
  EXPECT_TRUE(b[6]);
  EXPECT_FALSE(b[1]);
  EXPECT_FALSE(b[4]);
  EXPECT_FALSE(b[5]);
  EXPECT_FALSE(b[7]);
}

TEST(BitsetTest, SetResetFlip) {
  bs::bitset<8> b;
  b.set(2);
  EXPECT_TRUE(b[2]);
  b.set();
  for (int i = 0; i < 8; i++) {
    EXPECT_TRUE(b[i]);
  }
  b.reset(2);
  EXPECT_FALSE(b[2]);
  b.flip(3);
  EXPECT_FALSE(b[3]);
  b.set();
  b.flip();
  for (int i = 0; i < 8; i++) {
    EXPECT_FALSE(b[i]);
  }
}

TEST(BitsetTest, ToULong) {
  bs::bitset<8> b(0b1100ULL);
  unsigned long v = b.to_ulong();
  EXPECT_EQ(v, 0b1100UL);
}

TEST(BitsetTest, ToString) {
  bs::bitset<4> b(0b1010ULL);
  std::string s = b.to_string();
  EXPECT_EQ(s, "1010");
}

TEST(BitsetTest, Count) {
  bs::bitset<8> b(0b10101100ULL);
  EXPECT_EQ(b.count(), 4u);
}

TEST(BitsetTest, Size) {
  bs::bitset<8> b;
  EXPECT_EQ(b.size(), 8u);
}

TEST(BitsetTest, AnyAllNone) {
  bs::bitset<8> b;
  EXPECT_FALSE(b.any());
  EXPECT_FALSE(b.all());
  EXPECT_TRUE(b.none());
  b.set();
  EXPECT_TRUE(b.any());
  EXPECT_TRUE(b.all());
  EXPECT_FALSE(b.none());
  b.reset();
  b.set(3);
  EXPECT_TRUE(b.any());
  EXPECT_FALSE(b.all());
  EXPECT_FALSE(b.none());
}

TEST(BitsetTest, EqualityOperator) {
  bs::bitset<8> b1(0b10101010ULL);
  bs::bitset<8> b2(0b10101010ULL);
  bs::bitset<8> b3(0b11110000ULL);
  EXPECT_TRUE(b1 == b2);
  EXPECT_FALSE(b1 == b3);
}

TEST(BitsetTest, OperatorAND) {
  bs::bitset<8> b1(0b11110000ULL);
  bs::bitset<8> b2(0b10101010ULL);
  auto b3 = b1 & b2;
  EXPECT_EQ(b3.to_ulong(), 0b10100000UL);
}

TEST(BitsetTest, OperatorOR) {
  bs::bitset<8> b1(0b11110000ULL);
  bs::bitset<8> b2(0b10101010ULL);
  auto b3 = b1 | b2;
  EXPECT_EQ(b3.to_ulong(), 0b11111010UL);
}

TEST(BitsetTest, OperatorXOR) {
  bs::bitset<8> b1(0b11110000ULL);
  bs::bitset<8> b2(0b10101010ULL);
  auto b3 = b1 ^ b2;
  EXPECT_EQ(b3.to_ulong(), 0b01011010UL);
}

TEST(BitsetTest, ShiftLeft) {
  bs::bitset<8> b(0b10110000ULL);
  b <<= 2;
  EXPECT_EQ(b.to_ulong(), 0b11000000UL);
}

TEST(BitsetTest, ShiftRight) {
  bs::bitset<8> b(0b10110000ULL);
  b >>= 3;
  ASSERT_EQ(b[0], 0);
  ASSERT_EQ(b[1], 1);
  ASSERT_EQ(b[2], 1);
  ASSERT_EQ(b[3], 0);
  ASSERT_EQ(b[4], 1);
  ASSERT_EQ(b[5], 0);
  ASSERT_EQ(b[6], 0);
  ASSERT_EQ(b[7], 0);
  EXPECT_EQ(b.to_ulong(), 0b00010110UL);
}

TEST(BitsetTest, OutOfRange) {
  bs::bitset<8> b;
  EXPECT_THROW(b.set(8), std::out_of_range);
  EXPECT_THROW(b.test(8), std::out_of_range);
}

using b4_t = bs::bitset<4>;

TEST(BitsetTest, RefGet) {
  b4_t b{1};
  b4_t::reference ref = b[0];
  EXPECT_EQ(static_cast<bool>(ref), true);
}

TEST(BitsetTest, RefAssign) {
  b4_t b{1};
  b4_t::reference ref = b[0];
  ASSERT_EQ(static_cast<bool>(ref), true);
  ref = false;
  ASSERT_EQ(static_cast<bool>(ref), false);
  EXPECT_EQ(b[0], false);
}

TEST(BitsetTest, RefCopyCtor) {
  b4_t b{1};
  b4_t::reference ref = b[0];
  b4_t::reference ref2 = ref;
  ASSERT_TRUE(static_cast<bool>(ref));
  EXPECT_TRUE(static_cast<bool>(ref2));
  ref = false;
  EXPECT_EQ(static_cast<bool>(ref), static_cast<bool>(ref2));
  EXPECT_EQ(static_cast<bool>(ref), b[0]);
  EXPECT_FALSE(static_cast<bool>(ref));
  ref = true;
  EXPECT_EQ(static_cast<bool>(ref), static_cast<bool>(ref2));
  EXPECT_EQ(static_cast<bool>(ref), b[0]);
  EXPECT_TRUE(static_cast<bool>(ref));
  b[0] = false;
  EXPECT_EQ(static_cast<bool>(ref), static_cast<bool>(ref2));
  EXPECT_EQ(static_cast<bool>(ref), b[0]);
  EXPECT_FALSE(static_cast<bool>(ref));
}

TEST(BitsetTest, RefCopyAssign) {
  b4_t b{1};
  b4_t::reference ref = b[0];
  b4_t::reference ref2 = b[1];
  ref2 = ref;
  ASSERT_TRUE(static_cast<bool>(ref));
  EXPECT_TRUE(static_cast<bool>(ref2));
  ref = false;
  EXPECT_EQ(static_cast<bool>(ref), static_cast<bool>(ref2));
  EXPECT_EQ(static_cast<bool>(ref), b[0]);
  EXPECT_FALSE(static_cast<bool>(ref));
  ref = true;
  EXPECT_EQ(static_cast<bool>(ref), static_cast<bool>(ref2));
  EXPECT_EQ(static_cast<bool>(ref), b[0]);
  EXPECT_TRUE(static_cast<bool>(ref));
  b[0] = false;
  EXPECT_EQ(static_cast<bool>(ref), static_cast<bool>(ref2));
  EXPECT_EQ(static_cast<bool>(ref), b[0]);
  EXPECT_FALSE(static_cast<bool>(ref));
}


TEST(BitsetTest, RefNot) {
  b4_t b{1};
  b4_t::reference ref = b[0];
  EXPECT_FALSE(~ref);
}

TEST(BitsetTest, RefFlip) {
  b4_t b{1};
  b4_t::reference ref = b[0];
  EXPECT_FALSE(static_cast<bool>(ref.flip()));
  EXPECT_FALSE(static_cast<bool>(ref));
  EXPECT_FALSE(b[0]);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
