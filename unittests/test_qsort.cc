#include <gtest/gtest.h>
#include <gtest/test.hpp>

#ifdef _OPENMP
#include <mummer/openmp_qsort.hpp>

// #ifndef _OPENMP
// #error OPENMP is required
// #endif

namespace {
TEST(Qsort, Integers) {
  static size_t size = 100000;
  std::uniform_int_distribution<int> randnb(-1000000, 1000000);

  std::vector<int> numbers;
  for(size_t i = 0; i < size; ++i)
    numbers.push_back(randnb(rand_gen));

  EXPECT_FALSE(std::is_sorted(numbers.cbegin(), numbers.cend()));
  openmp_qsort(numbers.begin(), numbers.end());
  EXPECT_TRUE(std::is_sorted(numbers.cbegin(), numbers.cend()));
}
} // empty namespace

#endif
