#include <gtest/gtest.h>
#include <gtest/test.hpp>
#include <algorithm>

#include <mummer/sparseSA.hpp>

namespace {
TEST(SparseSA, ComputeLCP) {
  const std::string seq = sequence(1000) + std::string(300, 'N') + sequence(400) + std::string(400, 'N') + sequence(500);
  const long        N   = seq.size();
  const auto        sa  = mummer::mummer::sparseSA::create_auto(seq.c_str(), N, 10, true);

  EXPECT_EQ((unsigned int)0, sa.LCP[0]);
  for(long i = 1; i < N; ++i) {
    SCOPED_TRACE(::testing::Message() << "i:" << i);
    EXPECT_GE(N, sa.SA[i] + sa.LCP[i]);
    EXPECT_GE(N, sa.SA[i - 1] + sa.LCP[i]);
    EXPECT_LE((unsigned int)0, sa.LCP[i]);
    EXPECT_EQ(seq.substr(sa.SA[i-1], sa.LCP[i]), seq.substr(sa.SA[i], sa.LCP[i]));
    SCOPED_TRACE(::testing::Message() << "SA[i-1]:" << sa.SA[i-1] << " SA[i]:" << sa.SA[i] << " LCP[i]:" << sa.LCP[i]);
    if(sa.SA[i - 1] + sa.LCP[i] + 1 < N && sa.SA[i] + sa.LCP[i] + 1 < N) {
      EXPECT_NE(seq[sa.SA[i - 1] + sa.LCP[i]], seq[sa.SA[i] + sa.LCP[i]]);
    }
  }
}

void compareSA(const mummer::mummer::sparseSA& sa, const mummer::mummer::sparseSA& sa2) {
  EXPECT_EQ(sa._4column, sa2._4column);
  EXPECT_EQ(sa.K, sa2.K);
  EXPECT_EQ(sa.S.al_, sa2.S.al_);
  EXPECT_EQ(sa.S.l_, sa2.S.l_);
  EXPECT_EQ(0, memcmp(sa.S.s_, sa2.S.s_, sa.S.al_));
  EXPECT_EQ(sa.N, sa2.N);
  EXPECT_EQ(sa.logN, sa2.logN);
  EXPECT_EQ(sa.NKm1, sa2.NKm1);
  EXPECT_EQ(sa.SA.is_small, sa2.SA.is_small);
  EXPECT_EQ(sa.SA.size(), sa2.SA.size());
  if(sa.SA.is_small)
    EXPECT_TRUE(std::equal(sa.SA.small.cbegin(), sa.SA.small.cend(), sa2.SA.small.cbegin()));
  else
    EXPECT_TRUE(std::equal(sa.SA.large.cbegin(), sa.SA.large.cend(), sa2.SA.large.cbegin()));
  EXPECT_EQ(sa.ISA.is_small, sa2.SA.is_small);
  EXPECT_EQ(sa.ISA.size(), sa2.SA.size());
  if(sa.ISA.is_small)
    EXPECT_TRUE(std::equal(sa.ISA.small.cbegin(), sa.ISA.small.cend(), sa2.ISA.small.cbegin()));
  else
    EXPECT_TRUE(std::equal(sa.ISA.large.cbegin(), sa.ISA.large.cend(), sa2.ISA.large.cbegin()));
  EXPECT_TRUE(std::equal(sa.LCP.vec.cbegin(), sa.LCP.vec.cend(), sa2.LCP.vec.cbegin()));
  EXPECT_TRUE(std::equal(sa.LCP.M.cbegin(), sa.LCP.M.cend(), sa2.LCP.M.cend()));
  EXPECT_EQ(&sa2.SA, sa2.LCP.sa);
  EXPECT_TRUE(std::equal(sa.CHILD.cbegin(), sa.CHILD.cend(), sa2.CHILD.cbegin()));
  EXPECT_TRUE(std::equal(sa.KMR.cbegin(), sa.KMR.cend(), sa2.KMR.cbegin()));
  EXPECT_EQ(sa.hasChild, sa2.hasChild);
  EXPECT_EQ(sa.hasSufLink, sa2.hasSufLink);
  EXPECT_EQ(sa.hasKmer, sa2.hasKmer);
  EXPECT_EQ(sa.kMerSize, sa2.kMerSize);
  EXPECT_EQ(sa.kMerTableSize, sa2.kMerTableSize);
  EXPECT_EQ(sa.sparseMult, sa2.sparseMult);
  EXPECT_EQ(sa.nucleotidesOnly, sa2.nucleotidesOnly);
}

class SparseSATest : public ::testing::TestWithParam<bool> { };

TEST_P(SparseSATest, SaveLoad) {
  SCOPED_TRACE(::testing::Message() << (GetParam() ? "Large" : "Small") << " SA");
  const std::string seq = sequence(10000);

  prefix_unlink prefix("test_save");

  const auto sa = mummer::mummer::sparseSA::create_auto(seq.c_str(), seq.size(), 10, true, 1, GetParam());
  ASSERT_TRUE(sa.save(prefix.path));
  mummer::mummer::sparseSA sa2(seq.c_str(), seq.size(), prefix.path);
  { SCOPED_TRACE(::testing::Message() << "Loaded SA");
    compareSA(sa, sa2);
  }

  const auto sa3 = std::move(sa2);
  { SCOPED_TRACE(::testing::Message() << "Moved SA");
    compareSA(sa, sa3);
  }
} // SparseSA.SaveLoad
INSTANTIATE_TEST_CASE_P(SparseSA, SparseSATest, ::testing::Bool());
} // empty namespace
