#include <random>
#include <fstream>
#include <gtest/gtest.h>
#include <gtest/test.hpp>
#include <mummer/nucmer.hpp>

namespace {
// char comp(char b) {
//   switch(b) {
//   case 'a': return 't'; case 'A': return 'T';
//   case 'c': return 'g'; case 'C': return 'G';
//   case 'g': return 'c'; case 'G': return 'C';
//   case 't': return 'a'; case 'T': return 'A';
//   default: return 'n';
//   }
// }

// void assert_good_alignment(const mummer::postnuc::Alignment& al, const std::string& s1, const std::string& s2) {
//   SCOPED_TRACE(::testing::Message() << "al:" << al);
//   unsigned int errors = 0, indels = 0;
//   const auto it_end = mummer::postnuc::error_iterator_type(al, s1.c_str());
//   auto it = mummer::postnuc::error_iterator_type(al, s1.c_str(), s2.c_str(), s2.length());

//   EXPECT_LE(al.sA - 1, it->posA);
//   ASSERT_LE(0, it->posA);
//   ASSERT_GT(s1.size(), (size_t)it->posA);
//   ASSERT_LE(0, it->posB);
//   ASSERT_GT(s2.size(), (size_t)it->posB);

//   for( ; it != it_end; ++it) {
//     SCOPED_TRACE(::testing::Message() << "it<type:" << it->type << ",posA:" << it->posA << ",posB:" << it->posB << '>');
//     EXPECT_NE(mummer::postnuc::NONE, it->type);
//     ++errors;
//     indels += it->type != mummer::postnuc::MISMATCH;
//     ASSERT_LE(0, it->posA);
//     ASSERT_GT(s1.size(), (size_t)it->posA);
//     ASSERT_LE(0, it->posB);
//     ASSERT_GT(s2.size(), (size_t)it->posB);
//     if(it->type == mummer::postnuc::MISMATCH) {
//       EXPECT_NE(it->baseA, it->baseB);
//       EXPECT_EQ(s1[it->posA], it->baseA);
//       EXPECT_EQ(s2[it->posB], al.dirB == 1 ? it->baseB : comp(it->baseB));
//     }
//   }
// //   EXPECT_EQ(al.eA + 1, i);
// //   EXPECT_EQ(al.dirB == 1 ? al.eB + 1: (s2.size() - al.eB), j);
//   EXPECT_EQ(al.SimErrors, errors);
//   EXPECT_EQ(al.delta.size(), indels);
//   EXPECT_EQ(al.eA - 1, it->posA);
// }

TEST(Nucmer, PairSequences) {
  std::string s1 = sequence(1000);
  std::string s2 = s1.substr(900) + sequence(900);
  s1.erase(950, 1); // indel
  s2.erase(75, 1);  // indel
  s2[25] = (s2[25] == 'a' ? 'c' : 'a'); // substitution
  EXPECT_EQ((size_t)999, s1.size());
  EXPECT_EQ((size_t)999, s2.size());
  std::istringstream refstream(s1);
  mummer::nucmer::Options opts;
  opts.minmatch(10).mincluster(15);

  const auto a = mummer::nucmer::align_sequences(s1.c_str(), s1.length(),
                                                 s2.c_str(), s2.length(), opts);

  EXPECT_GT(a.size(), (size_t)0);
  const auto al = std::find_if(a.cbegin(), a.cend(), [](const mummer::postnuc::Alignment& al) { return al.sA == 901 && al.eA == 999; });
  ASSERT_NE(a.cend(), al);
  EXPECT_EQ(901, al->sA);
  EXPECT_EQ(999, al->eA);
  EXPECT_EQ(1, al->sB);
  EXPECT_EQ(99, al->eB);
  EXPECT_EQ(1, al->dirB);
  EXPECT_EQ(3, al->Errors);
  EXPECT_EQ(3, al->SimErrors);
  EXPECT_EQ((size_t)2, al->delta.size());
  // for(const auto& al : a)
  //   assert_good_alignment(al, s1, s2);
}

TEST(Nucmer, LongSequences) {
  const std::string s1 = sequence(1000);
  const std::string s2 = s1.substr(900) + sequence(900);

  std::istringstream refstream(std::string(">ref\n") + s1);
  mummer::nucmer::Options opts;
  mummer::nucmer::FileAligner falign(refstream, opts);
  const mummer::nucmer::FastaRecordSeq query_record(s2, "query");

  //  std::vector<mummer::postnuc::Alignment> a2;
  int nb = 0;
  auto append_matches = [&](std::vector<mummer::postnuc::Alignment>&& al,
                       const mummer::nucmer::FastaRecordPtr& ref,
                       const mummer::nucmer::FastaRecordSeq& query) {
    ++nb;
  };
  falign.align_long_sequences(query_record, append_matches);
  EXPECT_GT(nb, 0);

} // Nucmer.LongSequences

} // empty namespace
