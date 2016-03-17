#include <thread>
#include <vector>

#include <gtest/gtest.h>
#include <jellyfish/cooperative_pool2.hpp>

namespace {
// Generate all numbers in [0, producers * max)
class sequence : public jellyfish::cooperative_pool2<sequence, int> {
  typedef jellyfish::cooperative_pool2<sequence, int> super;

  const uint32_t        max_;
  std::vector<uint32_t> cur_;

public:
  sequence(uint32_t producers, uint32_t threads, uint32_t max)
    : super(producers, 4 * threads)
    , max_(max)
    , cur_(producers, 0)
  { }

  bool produce(uint32_t i, int& e) {
    auto& cur = cur_[i];
    if(cur < max_) {
      e = i * max_ + cur++;
      return false;
    }
    return true;
  }
};

class CooperativePoolTest : public ::testing::TestWithParam<uint32_t> {
public:
  static const uint32_t max_ = 10000;
  std::vector<unsigned char> check_;
  unsigned int nb_threads_;
  sequence seq_;  // Generator of ints

  CooperativePoolTest()
    : check_(GetParam() * max_, 0)
    , nb_threads_(std::thread::hardware_concurrency())
    , seq_(GetParam(), nb_threads_, max_)

  { }

  bool check() const {
    for(const auto it : check_)
      if(it != 1)
        return false;
    return true;
  }

  bool check_print() const {
    if(check())
      return true;

    for(int it : check_)
      std::cout << it;
    std::cout << "\n----------\n";
    // for(auto it = seq_.check_.cbegin(); it != seq_.check_.cend(); ++it)
    //   std::cout << (int)*it;
    std::cout << "\n";
    return false;
  }

  //protected:
};

void thread_work(sequence* const seq, std::vector<unsigned char>* const check) {
  while(true) {
    sequence::job j(*seq);
    if(j.is_empty())
      break;
    ++(*check)[*j];
  }
}

TEST_P(CooperativePoolTest, Ints) {
  std::vector<std::thread> threads; // Thread handles

  for(unsigned i = 0; i < nb_threads_; ++i)
    threads.push_back(std::thread(thread_work, &seq_, &check_));
  for(auto& th : threads)
    th.join();
  EXPECT_TRUE(check_print());
}

INSTANTIATE_TEST_CASE_P(CooperativePool,
                        CooperativePoolTest,
                        ::testing::Range((uint32_t)1, (uint32_t)5));

} // namespace {
