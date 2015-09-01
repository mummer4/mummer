#include <random>
#include <fstream>
#include <gtest/gtest.h>
#include <mummer/nucmer.hpp>

namespace {
std::string sequence(size_t len) {
  static char bases[4] = { 'a', 'c', 'g', 't' };
  std::default_random_engine         gen;
  std::uniform_int_distribution<int> dist(0, 3);
  std::string                        result;
  for(size_t i = 0; i < len; ++i)
    result += bases[dist(gen)];
  return result;
}

// char comp(char b) {
//   switch(b) {
//   case 'a': case 'A': return 't';
//   case 'c': case 'C': return 'g';
//   case 'g': case 'G': return 'c';
//   case 't': case 'T': return 'a';
//   default: return 'n';
//   }
// }
// void assert_good_alignment(const mummer::postnuc::Alignment& al, const std::string& s1, const std::string& s2) {
//   unsigned int errors = 0;
//   auto i = al.sA;
//   auto j = (al.dirB == 1 ? al.sB : s2.size() - al.sB);
//   auto ni = i - 1;
//   if(!al.delta.empty())
//     ni += std::abs(al.delta[0]) - 1;
//   size_t k = 0;

//   while(i <= al.eA) {
//     const char b2 = (al.dirB == 1 ? s2[j - 1] : comp(s2[j - 1]));
//     if(i == ni && k < al.delta.size()) {
//       while(i == ni && k < al.delta.size()) {
//         if(al.delta[k] > 0)
//           i += 1;
//         else
//           j += al.dirB;
//         errors +=1;
//         k += 1;
//         ni = i;
//         if(k < al.delta.size())
//           ni += std::abs(al.delta[k]);
//       }
//     } else {
//       if(s1[i - 1] != b2)
//         errors += 1;
//       i += 1;
//       j += al.dirB;
//     }
//   }

//   EXPECT_EQ(al.eA + 1, i);
//   EXPECT_EQ(al.dirB == 1 ? al.eB + 1: (s2.size() - al.eB), j);
//   EXPECT_EQ(al.delta.size(), k);
//   EXPECT_EQ(al.SimErrors, errors);
// }

TEST(Nucmer, PairSequences) {
  std::string s1 = sequence(1000);
  std::string s2 = s1.substr(900) + sequence(900);
  s1.erase(950, 1); // indel
  s2.erase(75, 1);  // indel
  s2[25] = (s2[25] == 'A' ? 'C' : 'A'); // substitution
  EXPECT_EQ((size_t)999, s1.size());
  EXPECT_EQ((size_t)999, s2.size());

  { std::ofstream os("test1.fa");
    os << ">ref\n" << s1 << '\n';
  }
  { std::ofstream os("test2.fa");
    os << ">qry\n" << s2 << '\n';
  }
  mummer::nucmer::Options opts;
  opts.minmatch(10).mincluster(10);
  const auto a = mummer::nucmer::align_sequences(s1.c_str(), s1.length(),
                                                 s2.c_str(), s2.length(), opts);

  EXPECT_GT(a.size(), (size_t)0);
  const auto al = std::find_if(a.cbegin(), a.cend(), [](const mummer::postnuc::Alignment& al) { return al.sA == 901 && al.eA == 999; });
  EXPECT_NE(a.cend(), al);
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
} // empty namespace
