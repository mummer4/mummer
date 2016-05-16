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


#ifndef __SKIP_LIST_COMMON_HPP__
#define __SKIP_LIST_COMMON_HPP__

namespace mt_skip_list {
namespace imp {

// XOR RNG by George Marsalia
struct xor_random {
  typedef uint64_t rand_type;
  uint64_t x;
  // xor64
  static rand_type gen(rand_type& y) {
    y ^= y << 13;
    y ^= y >> 7;
    y ^= y << 17;
    return y;
  }
  rand_type operator()() {
    return gen(x);
  }

  // xor128
  // uint64_t x, y, z, w;
  // rand_type operator()() {
  //   uint64_t t = x ^ (x << 5);
  //   x = y;
  //   y = z;
  //   z = w;
  //   w = (w ^ (w >> 29)) ^ (t ^ (t >> 12));
  //   return w;
  // }
  explicit xor_random() : x(88172645463325252LL) { };
  explicit xor_random(uint64_t seed, int n = 10) : x(seed) {
    for(int i = 0; i < n; ++i)
      this->operator()();
  }
};

inline int ctz(unsigned int x) { return __builtin_ctz(x); }
inline int ctz(unsigned long x) { return __builtin_ctzl(x); }
inline int ctz(unsigned long long x) { return __builtin_ctzll(x); }

// Return the height of a tower of pointer. Specialized for 2 and 4.
template<typename Random, int p>
struct random_height;
template<typename Random>
struct random_height<Random, 2> {
  Random rng;
  static inline int gen(typename Random::rand_type x) {
    return (x == 0 ? 8*sizeof(typename Random::rand_type) : ctz(x)) + 1;
  }
  inline int operator()() {
    return gen(rng());
  }
  random_height(const Random& rng_ = Random()) : rng(rng_) { }
};
template<typename Random>
struct random_height<Random, 4> {
  Random rng;
  static inline int gen(typename Random::rand_type x) {
    return (x == 0 ? 4 * sizeof(typename Random::rand_type) : ctz(x) >> 1) + 1;
  }
  inline int operator()() {
    return gen(rng());
  }
  random_height(const Random& rng_ = Random()) : rng(rng_) { }
};

// Comparator for the first element of a pair
template<typename Pair, class Compare>
struct first_comp {
  typedef Pair pair_type;
  typedef Compare comp_type;
  Compare comp;
  bool operator()(const pair_type& p1, const pair_type& p2) const {
    return comp(p1.first, p2.first);
  }
  bool operator()(const typename pair_type::first_type& p1, const pair_type& p2) const {
    return comp(p1, p2.first);
  }
  bool operator()(const pair_type& p1, const typename pair_type::first_type& p2) const {
    return comp(p1.first, p2);
  }
  first_comp(const Compare& comp_ = Compare()) : comp(comp_) { }
};

// Upper bound on height, depending on p.
template<int p>
struct height_bound { };
template<>
struct height_bound<2> { static const int value = 64; };
template<>
struct height_bound<4> { static const int value = 32; };

template<typename T>
struct block_allocator {
  struct block {
    block*              next;
    std::atomic<size_t> used;
    block() : next(nullptr), used(0) { }
  };

  block* head;
};
} // namespace mt_skip_list
} // namespace imp

#endif /* __SKIP_LIST_COMMON_HPP__ */
