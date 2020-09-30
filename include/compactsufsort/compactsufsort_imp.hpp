#ifndef __COMPACTSURFSORT_IMP_H__
#define __COMPACTSURFSORT_IMP_H__

#include <iostream>
#include <memory>
#include "divsufsort_private.h"
#include "sssort_imp.hpp"
#include "trsort_imp.hpp"
#include <mummer/timer.hpp>

#undef _OPENMP
#ifdef _OPENMP
# include <omp.h>
#endif


namespace compactsufsort_imp {
template<typename CHARPTR, typename SAIDPTR>
struct SA {
  typedef typename type_traits<SAIDPTR>::SAIDX          SAIDX;
  typedef typename const_iterator_traits<SAIDPTR>::type CSAIDPTR;
  static const size_t                                   ALPHABET_SIZE = alphabet_traits<CHARPTR>::size;

  inline static SAIDX& bucket(SAIDX* b, saint_t c0) {
    return b[c0];
  }
  inline static SAIDX& bucket(SAIDX* b, saint_t c0, saint_t c1) {
    return b[c1 * ALPHABET_SIZE + c0];
  }
  inline static SAIDX& bucket_star(SAIDX* b, saint_t c0, saint_t c1) {
    return b[c0 * ALPHABET_SIZE + c1];
  }


  /* Sorts suffixes of type B*. */
  static SAIDX
  sort_typeBstar(CHARPTR T, SAIDPTR SA,
                 SAIDX *bucket_A, SAIDX *bucket_B,
                 SAIDX n) {
    SAIDPTR PAb, ISAb, buf;

#ifdef _OPENMP
    SAIDPTR curbuf;
    SAIDX l;
#endif
    SAIDX i, j, k, t, m, bufsize;
    saint_t c0, c1;
#ifdef _OPENMP
    saint_t d0, d1;
    int tmp;
#endif

    /* Initialize bucket arrays. */
    for(SAIDX i = 0; i < (SAIDX)ALPHABET_SIZE; ++i) { bucket_A[i] = 0; }
    for(SAIDX i = 0; i < (SAIDX)ALPHABET_SIZE * (SAIDX)ALPHABET_SIZE; ++i) { bucket_B[i] = 0; }

    /* Count the number of occurrences of the first one or two characters of each
       type A, B and B* suffix. Moreover, store the beginning position of all
       type B* suffixes into the array SA. */
    { TIME_SCOPE("Count suffixes");
    for(i = n - 1, m = n, c0 = T[n - 1]; 0 <= i;) {
      /* type A suffix. */
      do { ++bucket(bucket_A, c1 = c0); } while((0 <= --i) && ((c0 = T[i]) >= c1));
      if(0 <= i) {
        /* type B* suffix. */
        ++bucket_star(bucket_B, c0, c1);
        SA[--m] = i;
        /* type B suffix. */
        for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1); --i, c1 = c0) {
          ++bucket(bucket_B, c0, c1);
        }
      }
    }
    m = n - m;
    } // time_scope
    /*
      note:
      A type B* suffix is lexicographically smaller than a type B suffix that
      begins with the same first two characters.
    */

    /* Calculate the index of start/end point of each bucket. */
    { TIME_SCOPE("calculate start/end");
    for(c0 = 0, i = 0, j = 0; c0 < (saint_t)ALPHABET_SIZE; ++c0) {
      t = i + bucket(bucket_A, c0);
      bucket(bucket_A, c0) = i + j; /* start point */
      i = t + bucket(bucket_B, c0, c0);
      for(c1 = c0 + 1; c1 < (saint_t)ALPHABET_SIZE; ++c1) {
        j += bucket_star(bucket_B, c0, c1);
        bucket_star(bucket_B, c0, c1) = j; /* end point */
        i += bucket(bucket_B, c0, c1);
      }
    }
    } // time_scope

    if(0 < m) {
      { TIME_SCOPE("sort 2 char");
      /* Sort the type B* suffixes by their first two characters. */
      PAb = SA + n - m; ISAb = SA + m;
      for(SAIDX i = m - 2; 0 <= i; --i) {
        t = PAb[i], c0 = T[t], c1 = T[t + 1];
        SA[--bucket_star(bucket_B, c0, c1)] = i;
      }
      t = PAb[m - 1], c0 = T[t], c1 = T[t + 1];
      SA[--bucket_star(bucket_B, c0, c1)] = m - 1;
      } // time_scope

      /* Sort the type B* substrings using sssort. */
      { TIME_SCOPE("sssort");
#ifdef _OPENMP
      tmp = omp_get_max_threads();
      buf = SA + m, bufsize = (n - (2 * m)) / tmp;
      c0 = ALPHABET_SIZE - 2, c1 = ALPHABET_SIZE - 1, j = m;
#pragma omp parallel default(shared) private(curbuf, k, l, d0, d1, tmp)
      {
        tmp = omp_get_thread_num();
        curbuf = buf + tmp * bufsize;
        k = 0;
        for(;;) {
#pragma omp critical(sssort_lock)
          {
            if(0 < (l = j)) {
              d0 = c0, d1 = c1;
              do {
                k = bucket_star(bucket_B, d0, d1);
                if(--d1 <= d0) {
                  d1 = ALPHABET_SIZE - 1;
                  if(--d0 < 0) { break; }
                }
              } while(((l - k) <= 1) && (0 < (l = k)));
              c0 = d0, c1 = d1, j = k;
            }
          }
          if(l == 0) { break; }
          ss<CHARPTR, SAIDPTR>::sort(T, PAb, SA + k, SA + l,
                            curbuf, bufsize, (SAIDX)2, n, *(SA + k) == (m - 1));
        }
      }
#else
      buf = SA + m, bufsize = n - (2 * m);
      for(c0 = ALPHABET_SIZE - 2, j = m; 0 < j; --c0) {
        for(c1 = ALPHABET_SIZE - 1; c0 < c1; j = i, --c1) {
          i = bucket_star(bucket_B, c0, c1);
          if(1 < (j - i)) {
            ss<CHARPTR, SAIDPTR>::sort(T, PAb, SA + i, SA + j,
                              buf, bufsize, (SAIDX)2, n, *(SA + i) == (m - 1));
          }
        }
      }
#endif
      } // time_scope

      /* Compute ranks of type B* substrings. */
      { TIME_SCOPE("Ranks B*");
      for(SAIDX i = m - 1; 0 <= i; --i) {
        if(0 <= SA[i]) {
          j = i;
          do { ISAb[SA[i]] = i; } while((0 <= --i) && (0 <= SA[i]));
          SA[i + 1] = i - j;
          if(i <= 0) { break; }
        }
        j = i;
        do { ISAb[SA[i] = ~SA[i]] = j; } while(SA[--i] < 0);
        ISAb[SA[i]] = j;
      }
      } // time_scope

      /* Construct the inverse suffix array of type B* suffixes using trsort. */
      { TIME_SCOPE("ISA");
      tr<SAIDPTR>::sort(ISAb, SA, m, (SAIDX)1);
      } // time_scope

      /* Set the sorted order of tyoe B* suffixes. */
      { TIME_SCOPE("sorted order B*");
      for(i = n - 1, j = m, c0 = T[n - 1]; 0 <= i;) {
        for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) >= c1); --i, c1 = c0) { }
        if(0 <= i) {
          t = i;
          for(--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1); --i, c1 = c0) { }
          SA[ISAb[--j]] = ((t == 0) || (1 < (t - i))) ? t : ~t;
        }
      }
      } // time_scope

      /* Calculate the index of start/end point of each bucket. */
      { TIME_SCOPE("start/end each bucket");
      bucket(bucket_B, ALPHABET_SIZE - 1, ALPHABET_SIZE - 1) = n; /* end point */
      for(c0 = ALPHABET_SIZE - 2, k = m - 1; 0 <= c0; --c0) {
        i = bucket(bucket_A, c0 + 1) - 1;
        for(c1 = ALPHABET_SIZE - 1; c0 < c1; --c1) {
          t = i - bucket(bucket_B, c0, c1);
          bucket(bucket_B, c0, c1) = i; /* end point */

          /* Move all type B* suffixes to the correct position. */
          for(i = t, j = bucket_star(bucket_B, c0, c1);
              j <= k;
              --i, --k) { SA[i] = SA[k]; }
        }
        bucket_star(bucket_B, c0, c0 + 1) = i - bucket(bucket_B, c0, c0) + 1; /* start point */
        bucket(bucket_B, c0, c0) = i; /* end point */
      }
      } // time_scope
    }

    return m;
  }

  /* Constructs the suffix array by using the sorted order of type B* suffixes. */
  static void
  construct_SA(CHARPTR T, SAIDPTR SA,
               SAIDX *bucket_A, SAIDX *bucket_B,
               SAIDX n, SAIDX m) {
    SAIDPTR i, j, k = nullptr;
    SAIDX s;
    saint_t c0, c1, c2;

    if(0 < m) {
      /* Construct the sorted order of type B suffixes by using
         the sorted order of type B* suffixes. */
      for(c1 = ALPHABET_SIZE - 2; 0 <= c1; --c1) {
        /* Scan the suffix array from right to left. */
        for(i = SA + bucket_star(bucket_B, c1, c1 + 1),
              j = SA + bucket(bucket_A, c1 + 1) - 1, c2 = -1;
            i <= j;
            --j) {
          if(0 < (s = *j)) {
            assert(T[s] == c1);
            assert(((s + 1) < n) && (T[s] <= T[s + 1]));
            assert(T[s - 1] <= T[s]);
            *j = ~s;
            c0 = T[--s];
            if((0 < s) && (T[s - 1] > c0)) { s = ~s; }
            if(c0 != c2) {
              if(0 <= c2) { bucket(bucket_B, c2, c1) = k - SA; }
              k = SA + bucket(bucket_B, c2 = c0, c1);
            }
            assert(k < j);
            *k-- = s;
          } else {
            assert(((s == 0) && (T[s] == c1)) || (s < 0));
            *j = ~s;
          }
        }
      }
    }

    /* Construct the suffix array by using
       the sorted order of type B suffixes. */
    k = SA + bucket(bucket_A, c2 = T[n - 1]);
    *k++ = (T[n - 2] < c2) ? ~(n - 1) : (n - 1);
    /* Scan the suffix array from left to right. */
    for(i = SA, j = SA + n; i < j; ++i) {
      if(0 < (s = *i)) {
        assert(T[s - 1] >= T[s]);
        c0 = T[--s];
        if((s == 0) || (T[s - 1] < c0)) { s = ~s; }
        if(c0 != c2) {
          bucket(bucket_A, c2) = k - SA;
          k = SA + bucket(bucket_A, c2 = c0);
        }
        assert(i < k);
        *k++ = s;
      } else {
        assert(s < 0);
        *i = ~s;
      }
    }
  }

  static saint_t
  create(CHARPTR T, SAIDPTR SA, SAIDX n) {
    std::unique_ptr<SAIDX[]> bucket_A, bucket_B;
    SAIDX m;
    saint_t err = 0;

    /* Check arguments. */
    if((T == nullptr) || (n < 0)) { return -1; }
    else if(n == 0) { return 0; }
    else if(n == 1) { SA[0] = 0; return 0; }
    else if(n == 2) { m = (T[0] < T[1]); SA[m ^ 1] = 0, SA[m] = 1; return 0; }

    bucket_A.reset(new(std::nothrow) SAIDX[ALPHABET_SIZE]);
    bucket_B.reset(new(std::nothrow) SAIDX[ALPHABET_SIZE * ALPHABET_SIZE]);

    /* Suffix sort. */
    if(bucket_A && bucket_B) {
      m = sort_typeBstar(T, SA, bucket_A.get(), bucket_B.get(), n);
      construct_SA(T, SA, bucket_A.get(), bucket_B.get(), n, m);
    } else {
      err = -2;
    }

    return err;
  }

  template<typename ITERATOR>
  static int
  _compare(ITERATOR T, SAIDX Tsize,
           ITERATOR P, SAIDX Psize,
           SAIDX suf, SAIDX *match) {
    SAIDX i, j;
    saint_t r;
    for(i = suf + *match, j = *match, r = 0;
        (i < Tsize) && (j < Psize) && ((r = T[i] - P[j]) == 0); ++i, ++j) { }
    *match = j;
    return (r == 0) ? -(j != Psize) : r;
  }

  /* Search for the pattern P in the string T. */
  static std::pair<SAIDX, SAIDX>
  search(CHARPTR T, SAIDX Tsize, CSAIDPTR SA, SAIDX SAsize,
         CHARPTR P, SAIDX Psize) {
    SAIDX size, lsize, rsize, half;
    SAIDX match, lmatch, rmatch;
    SAIDX llmatch, lrmatch, rlmatch, rrmatch;
    SAIDX i, j, k;
    saint_t r;

    if((T == nullptr) || (P == nullptr) || (SA == nullptr) ||
       (Tsize < 0) || (Psize < 0)) { return std::make_pair((SAIDX)-1, (SAIDX)-1); }
    if(Tsize == 0) { return std::make_pair((SAIDX)0, (SAIDX)-1); }
    if(Psize == 0) { return std::make_pair(Tsize, (SAIDX)0); }

    for(i = j = k = 0, lmatch = rmatch = 0, size = SAsize, half = size >> 1;
        0 < size;
        size = half, half >>= 1) {
      match = std::min(lmatch, rmatch);
      r = _compare(T, Tsize, P, Psize, (SAIDX)SA[i + half], &match);
      if(r < 0) {
        i += half + 1;
        half -= (size & 1) ^ 1;
        lmatch = match;
      } else if(r > 0) {
        rmatch = match;
      } else {
        lsize = half, j = i, rsize = size - half - 1, k = i + half + 1;

        /* left part */
        for(llmatch = lmatch, lrmatch = match, half = lsize >> 1;
            0 < lsize;
            lsize = half, half >>= 1) {
          lmatch = std::min(llmatch, lrmatch);
          r = _compare(T, Tsize, P, Psize, (SAIDX)SA[j + half], &lmatch);
          if(r < 0) {
            j += half + 1;
            half -= (lsize & 1) ^ 1;
            llmatch = lmatch;
          } else {
            lrmatch = lmatch;
          }
        }

        /* right part */
        for(rlmatch = match, rrmatch = rmatch, half = rsize >> 1;
            0 < rsize;
            rsize = half, half >>= 1) {
          rmatch = std::min(rlmatch, rrmatch);
          r = _compare(T, Tsize, P, Psize, (SAIDX)SA[k + half], &rmatch);
          if(r <= 0) {
            k += half + 1;
            half -= (rsize & 1) ^ 1;
            rlmatch = rmatch;
          } else {
            rrmatch = rmatch;
          }
        }

        break;
      }
    }

    return std::make_pair(k - j, 0 < (k - j) ? j : i);
  }

  /* Search for the character c in the string T. */
  static std::pair<SAIDX, SAIDX>
  search(CHARPTR T, SAIDX Tsize, CSAIDPTR SA, SAIDX SAsize,
         saint_t c) {
    SAIDX size, lsize, rsize, half;
    SAIDX i, j, k, p;
    saint_t r;

    //  if(idx != NULL) { *idx = -1; }
    if((T == NULL) || (SA == NULL) || (Tsize < 0)) { return std::make_pair((SAIDX)-1, (SAIDX)-1); }
    if(Tsize == 0) { return std::make_pair((SAIDX)0, (SAIDX)-1); }

    for(i = j = k = 0, size = SAsize, half = size >> 1;
        0 < size;
        size = half, half >>= 1) {
      p = SA[i + half];
      r = (p < Tsize) ? T[p] - c : -1;
      if(r < 0) {
        i += half + 1;
        half -= (size & 1) ^ 1;
      } else if(r == 0) {
        lsize = half, j = i, rsize = size - half - 1, k = i + half + 1;

        /* left part */
        for(half = lsize >> 1;
            0 < lsize;
            lsize = half, half >>= 1) {
          p = SA[j + half];
          r = (p < Tsize) ? T[p] - c : -1;
          if(r < 0) {
            j += half + 1;
            half -= (lsize & 1) ^ 1;
          }
        }

        /* right part */
        for(half = rsize >> 1;
            0 < rsize;
            rsize = half, half >>= 1) {
          p = SA[k + half];
          r = (p < Tsize) ? T[p] - c : -1;
          if(r <= 0) {
            k += half + 1;
            half -= (rsize & 1) ^ 1;
          }
        }

        break;
      }
    }

    return std::make_pair(k - j, 0 < (k - j) ? j : i);
  }

  static saint_t
  check(CHARPTR T, SAIDPTR SA,
        SAIDX n, saint_t verbose) {
    SAIDX C[ALPHABET_SIZE];
    SAIDX i, p, q, t;
    saint_t c;

    if(verbose) { std::cerr << "sufcheck: "; }

    /* Check arguments. */
    if((T == nullptr) || (n < 0)) {
      if(verbose) { std::cerr << "Invalid arguments.\n"; }
      return -1;
    }
    if(n == 0) {
      if(verbose) { std::cerr << "Done.\n"; }
      return 0;
    }

    /* check range: [0..n-1] */
    for(i = 0; i < n; ++i) {
      if((SA[i] < 0) || (n <= SA[i])) {
        if(verbose) {
          std::cerr << "Out of the range [0," << (n-1) <<  "].\n"
                    << "  SA[" << i << "]=" << SA[i] << "\n";
        }
        return -2;
      }
    }

    /* check first characters. */
    for(i = 1; i < n; ++i) {
      if(T[SA[i - 1]] > T[SA[i]]) {
        if(verbose) {
          std::cerr << "Suffixes in wrong order.\n"
                    << "  T[SA["  << (i-1) << "]=" << SA[i-1] << "]=" << T[SA[i-1]]
                    << " > T[SA[" << i     << "]=" << SA[i]   << "]=" << T[SA[i]]
                    << '\n';
        }
        return -3;
      }
    }

    /* check suffixes. */
    for(i = 0; i < (SAIDX)ALPHABET_SIZE; ++i) { C[i] = 0; }
    for(i = 0; i < n; ++i) { ++C[T[i]]; }
    for(i = 0, p = 0; i < (SAIDX)ALPHABET_SIZE; ++i) {
      t = C[i];
      C[i] = p;
      p += t;
    }

    q = C[T[n - 1]];
    C[T[n - 1]] += 1;
    for(i = 0; i < n; ++i) {
      p = SA[i];
      if(0 < p) {
        c = T[--p];
        t = C[c];
      } else {
        c = T[p = n - 1];
        t = q;
      }
      if((t < 0) || (p != SA[t])) {
        if(verbose) {
          std::cerr << "Suffix in wrong position.\n"
                    << "  SA[" << t << "]=" << ((0 <= t) ? SA[t] : -1) << " or\n"
                    << "  SA[" << i << "]=" << SA[i] << '\n';
        }
        return -4;
      }
      if(t != q) {
        ++C[c];
        if((n <= C[c]) || (T[SA[C[c]]] != c)) { C[c] = -1; }
      }
    }

    if(1 <= verbose) { std::cerr << "Done.\n"; }
    return 0;
  }

}; // struct SA

} // namespace compactsufsort_imp

#endif /* __COMPACTSURFSORT_IMP_H__ */
