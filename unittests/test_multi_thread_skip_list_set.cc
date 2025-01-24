/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


 #include <gtest/gtest.h>
#include <unittests/misc.hpp>
#include <mt_skip_list/set.hpp>

#include <string>
#include <vector>

namespace {
// Compute the power at compile time in a stupid way
template<long x, long n>
struct lpow {
  enum { value = x * lpow<x, n-1>::value };
};
template<long x>
struct lpow<x, 0> {
  enum { value = 1 };
};

typedef mt_skip_list::set<int> set_type;
typedef set_type::iterator set_iterator;
typedef set_type::const_iterator const_set_iterator;
typedef std::pair<set_iterator, bool> set_ins;

typedef std::set<int>::iterator std_iterator;
typedef std::set<int>::const_iterator const_std_iterator;
typedef std::pair<std_iterator, bool> std_ins;

TEST(MTSkipListSet, Init) {
  set_type sls;
  EXPECT_TRUE(sls.empty());
  EXPECT_EQ((size_t)0, sls.size());
  EXPECT_TRUE(sls.end() == sls.begin());
  //    EXPECT_EQ(((size_t)lpow<set_type::p, max_height>::value), sls.max_size());
}

TEST(MTSkipListSet, InitializerList) {
  const std::initializer_list<int> il = { 1, 3, -5, -7 };
  const set_type sls(il);
  std::set<int> set(il);
  EXPECT_EQ(il.size(), sls.size());
  EXPECT_TRUE(std::equal(set.cbegin(), set.cend(), sls.cbegin()));
} // MTSkipListSet.InitializerList

TEST(MTSkipListSet, Copy) {
  const std::initializer_list<int> il = { -1, 3, -5, 7, -9 };
  set_type sls(il);
  const std::set<int>set(il);

  set_type csls(sls);
  EXPECT_EQ(il.size(), csls.size());
  EXPECT_TRUE(std::equal(set.cbegin(), set.cend(), csls.cbegin()));

  set_type ccsls;
  ccsls = csls;
  EXPECT_EQ(il.size(), csls.size());
  EXPECT_TRUE(std::equal(set.cbegin(), set.cend(), ccsls.cbegin()));
} // MTSkipListSet.Copy

TEST(MTSkipListSet, Move) {
  const std::initializer_list<int> il = { -1, 3, -5, 7, -9 };
  set_type sls(il);
  const std::set<int>set(il);

  set_type csls(std::move(sls));
  EXPECT_EQ(il.size(), csls.size());
  EXPECT_TRUE(std::equal(set.cbegin(), set.cend(), csls.cbegin()));

  set_type ccsls;
  ccsls = std::move(csls);
  EXPECT_EQ(il.size(), ccsls.size());
  EXPECT_TRUE(std::equal(set.cbegin(), set.cend(), ccsls.cbegin()));
} // MTSkipListSet.Move


TEST(MTSkipListSet, InsertOneThread) {
  set_type         sls, slse;
  std::set<int>    set;
  EXPECT_TRUE(sls.empty());
  EXPECT_TRUE(set.empty());

  for(int i = 0; i < 100; ++i) {
    int     x        = random() % 50;
    std_ins set_res  = set.insert(x);
    set_ins sls_res  = sls.insert(x);
    set_ins sls_rese = slse.emplace(x);
    EXPECT_EQ(set_res.second, sls_res.second);
    EXPECT_EQ(set_res.second, sls_rese.second);
    EXPECT_EQ(*set_res.first, *sls_res.first);
    EXPECT_EQ(*set_res.first, *sls_rese.first);
    EXPECT_EQ(x, *sls_res.first);
    EXPECT_EQ(x, *sls_rese.first);
    EXPECT_EQ(set.empty(), sls.empty());
    EXPECT_EQ(set.empty(), slse.empty());
    EXPECT_EQ(set.size(), sls.size());
    EXPECT_EQ(set.size(), slse.size());
  }

  const_set_iterator set_it = sls.begin();
  const_set_iterator set_ite = slse.begin();
  for(const_std_iterator std_it = set.begin(); std_it != set.end();
      ++std_it, ++set_it, ++set_ite) {
    ASSERT_FALSE(set_it == sls.end());
    ASSERT_FALSE(set_ite == slse.end());
    EXPECT_EQ(*std_it, *set_it);
    EXPECT_EQ(*std_it, *set_ite);
  }
  ASSERT_TRUE(set_it == sls.end());
  ASSERT_TRUE(set_it == slse.end());
}

struct thread_insert_data {
  set_type         set;       // The set to add to
  std::atomic<int> ids;       // Counter to get a thread id
  std::vector<int> v;         // The data to insert
  std::atomic<int> new_elt;   // Nb new element inserted
  std::atomic<int> exist_elt; // Nb existing element inserted

  static const int         per_th = 10000; // Nb elements per thread
  thread_insert_data(int nb_threads)
    : set() //, std::less<int>(), xor_random(random()))
    , ids(-1)
    , new_elt(0), exist_elt(0)
  {
    const int n = per_th * nb_threads;
    for(int i = 0; i < n; ++i)
      v.push_back(random() % n);
  }
};
void* insert_from_thread(thread_insert_data* data) {
  auto& set = data->set;

  int tid = (data->ids += 1);
  int new_elt = 0, exist_elt = 0;
  for(int i = tid * data->per_th; i < (tid+1) * data->per_th; ++i) {
    auto res = set.insert(data->v[i]);
    if(res.second)
      ++new_elt;
    else
      ++exist_elt;
  }
  new_elt = (data->new_elt   += new_elt);
  exist_elt = (data->exist_elt += exist_elt);

  return 0;
}
TEST(MTSkipListSet, InsertManyThreads) {
  const int          nb_threads = 5;
  thread_insert_data data(nb_threads);

  pdo(nb_threads, insert_from_thread, &data);
  EXPECT_EQ(data.v.size(), (size_t)(data.new_elt + data.exist_elt));

  // Do the same single threads into a set and check the statistics
  std::set<int> std_set;
  int new_elt = 0, exist_elt = 0;
  for(auto it = data.v.begin(); it != data.v.end(); ++it) {
    auto res = std_set.insert(*it);
    if(res.second)
      ++new_elt;
    else
      ++exist_elt;
  }
  EXPECT_EQ(new_elt, data.new_elt);
  EXPECT_EQ(exist_elt, data.exist_elt);
  EXPECT_EQ(std_set.size(), data.set.size());

  for(int i = -1; i <= (int)data.v.size(); ++i) {
    EXPECT_EQ(std_set.count(i), data.set.count(i));
    { // Test find
      auto std_res = std_set.find(i);
      auto set_res = data.set.find(i);
      if(std_res == std_set.end())
        EXPECT_EQ(data.set.end(), set_res);
      else
        EXPECT_EQ(*std_res, *set_res);
    }
    { // Test equal_range
      auto std_res = std_set.equal_range(i);
      auto set_res = data.set.equal_range(i);
      EXPECT_EQ(std_res.first == std_res.second, set_res.first == set_res.second);
      EXPECT_EQ(std_res.second == std_set.end(), set_res.second == data.set.end());
      if(std_res.first != std_res.second) {
        EXPECT_EQ(i, *set_res.first);
        EXPECT_EQ(set_res.second, ++set_res.first);
      }
    }
    { // Test lower_bound
      auto std_res = std_set.lower_bound(i);
      auto set_res = data.set.lower_bound(i);
      if(std_res == std_set.end())
        EXPECT_EQ(data.set.end(), set_res);
      else
        EXPECT_EQ(*std_res, *set_res);
    }
    { // Test upper_bound
      auto std_res = std_set.upper_bound(i);
      auto set_res = data.set.upper_bound(i);
      if(std_res == std_set.end())
        EXPECT_EQ(data.set.end(), set_res);
      else
        EXPECT_EQ(*std_res, *set_res);
    }
  }
}

// Insert many overlapping value. Should always only be 1 in set.
void insert_overlapping(mt_skip_list::set<long>* set, int max, std::atomic<int>* ready) {
  --*ready;
  while(ready->load())
    //    std::this_thread::yield();
    ;
  for(int i = 0; i < max; ++i)
    set->insert(i);
}


TEST(MTSkipListSet, InsertOverlapping) {
  const int nb_threads = 5;
  const int max        = 10000000;
  mt_skip_list::set<long> set;
  std::atomic<int>        ready(nb_threads);

  { std::vector<std::thread> threads;
    for(int i = 0; i < nb_threads; ++i)
      threads.push_back(std::thread(insert_overlapping, &set, max, &ready));
    for(auto& th : threads)
      th.join();
  }

  EXPECT_EQ((size_t)max, set.size());
  auto it = set.cbegin();
  for(int i = 0; i < max; ++i, ++it) {
    EXPECT_FALSE(set.cend() == it);
    EXPECT_EQ(i, *it);
  }
}

} // namespace
