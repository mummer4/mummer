#ifndef __SSSORT_IMP_H__
#define __SSSORT_IMP_H__

#include <utility>
#include <algorithm>
using std::swap;

#include "divsufsort_private.h"
#include "const_iterator_traits.hpp"

namespace compactsufsort_imp {
  /* Compares two suffixes. */
template<typename CHARPTR, typename SAIDPTR1, typename SAIDPTR2>
typename std::enable_if<std::is_same<const char*, CHARPTR>::value, saint_t>::type
compare(CHARPTR T, SAIDPTR1 p1, SAIDPTR2 p2, typename type_traits<SAIDPTR1>::SAIDX depth) {
  auto l = std::min(*(p1 + 1), *(p2 + 1)) + 2;
  return memcmp(T + depth + *p1, T + depth + *p2, l > depth ? l - depth : 0);
}

template<typename CHARPTR, typename SAIDPTR1, typename SAIDPTR2>
typename std::enable_if<!std::is_same<const char*, CHARPTR>::value, saint_t>::type
compare(CHARPTR T, SAIDPTR1 p1, SAIDPTR2 p2, typename type_traits<SAIDPTR1>::SAIDX depth) {
  CHARPTR U1, U2, U1n, U2n;

  for(U1 = T + depth + *p1,
        U2 = T + depth + *p2,
        U1n = T + *(p1 + 1) + 2,
        U2n = T + *(p2 + 1) + 2;
      (U1 < U1n) && (U2 < U2n) && (*U1 == *U2);
      ++U1, ++U2) {
  }

  return U1 < U1n ? (U2 < U2n ? *U1 - *U2 : 1) : (U2 < U2n ? -1 : 0);
}

template<typename CHARPTR, typename SAIDPTR>
struct ss {
  typedef typename type_traits<SAIDPTR>::SAIDX SAIDX;
  typedef typename const_iterator_traits<SAIDPTR>::type CSAIDPTR;

#if (SS_BLOCKSIZE != 1) && (SS_INSERTIONSORT_THRESHOLD != 1)

  /* Insertionsort for small size groups */

  static void
  insertionsort(CHARPTR T, CSAIDPTR PA,
                   SAIDPTR first, SAIDPTR last, SAIDX depth) {
    SAIDPTR i, j;
    SAIDX t;
    saint_t r;

    for(i = last - 2; first <= i; --i) {
      for(t = *i, j = i + 1; 0 < (r = compare(T, PA + t, PA + *j, depth));) {
        do { *(j - 1) = *j; } while((++j < last) && (*j < 0));
        if(last <= j) { break; }
      }
      if(r == 0) { *j = ~*j; }
      *(j - 1) = t;
    }
  }

#endif /* (SS_BLOCKSIZE != 1) && (SS_INSERTIONSORT_THRESHOLD != 1) */


  /*---------------------------------------------------------------------------*/

#if (SS_BLOCKSIZE == 0) || (SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE)

  static void
  fixdown(CHARPTR Td, CSAIDPTR PA,
             SAIDPTR SA, SAIDX i, SAIDX size) {
    SAIDX j, k;
    SAIDX v;
    saint_t c, d, e;

    for(v = SA[i], c = Td[PA[v]]; (j = 2 * i + 1) < size; SA[i] = SA[k], i = k) {
      d = Td[PA[SA[k = j++]]];
      if(d < (e = Td[PA[SA[j]]])) { k = j; d = e; }
      if(d <= c) { break; }
    }
    SA[i] = v;
  }

  /* Simple top-down heapsort. */
  static void
  heapsort(CHARPTR Td, CSAIDPTR PA, SAIDPTR SA, SAIDX size) {
    SAIDX i, m;
    SAIDX t;

    m = size;
    if((size % 2) == 0) {
      m--;
      if(Td[PA[SA[m / 2]]] < Td[PA[SA[m]]]) { swap(SA[m], SA[m / 2]); }
    }

    for(i = m / 2 - 1; 0 <= i; --i) { fixdown(Td, PA, SA, i, m); }
    if((size % 2) == 0) { swap(SA[0], SA[m]); fixdown(Td, PA, SA, (SAIDX)0, m); }
    for(i = m - 1; 0 < i; --i) {
      t = SA[0], SA[0] = SA[i];
      fixdown(Td, PA, SA, (SAIDX)0, i);
      SA[i] = t;
    }
  }


  /*---------------------------------------------------------------------------*/

  /* Returns the median of three elements. */
  static SAIDPTR
  median3(CHARPTR Td, CSAIDPTR PA,
          SAIDPTR v1, SAIDPTR v2, SAIDPTR v3) {
    const auto x1 = Td[PA[*v1]], x2 = Td[PA[*v2]], x3 = Td[PA[*v3]];
    if(x1 < x2)
      return x2 < x3 ? v2 : v3;
    return x1 < x3 ? v1 : v3;
  }

  /* Returns the median of five elements. */
  static SAIDPTR
  median5(CHARPTR Td, CSAIDPTR PA,
          SAIDPTR v1, SAIDPTR v2, SAIDPTR v3, SAIDPTR v4, SAIDPTR v5) {
    if(Td[PA[*v2]] > Td[PA[*v3]]) { swap(v2, v3); }
    if(Td[PA[*v4]] > Td[PA[*v5]]) { swap(v4, v5); }
    if(Td[PA[*v2]] > Td[PA[*v4]]) { swap(v2, v4); swap(v3, v5); }
    if(Td[PA[*v1]] > Td[PA[*v3]]) { swap(v1, v3); }
    if(Td[PA[*v1]] > Td[PA[*v4]]) { swap(v1, v4); swap(v3, v5); }
    if(Td[PA[*v3]] > Td[PA[*v4]]) { return v4; }
    return v3;
  }

  /* Returns the pivot element. */
  static SAIDPTR
  pivot(CHARPTR Td, CSAIDPTR PA, SAIDPTR first, SAIDPTR last) {
    SAIDPTR middle;
    SAIDX t;

    t = last - first;
    middle = first + t / 2;

    if(t <= 512) {
      if(t <= 32) {
        return median3(Td, PA, first, middle, last - 1);
      } else {
        t >>= 2;
        return median5(Td, PA, first, first + t, middle, last - 1 - t, last - 1);
      }
    }
    t >>= 3;
    first  = median3(Td, PA, first, first + t, first + (t << 1));
    middle = median3(Td, PA, middle - t, middle, middle + t);
    last   = median3(Td, PA, last - 1 - (t << 1), last - 1 - t, last - 1);
    return median3(Td, PA, first, middle, last);
  }


  /*---------------------------------------------------------------------------*/

  /* Binary partition for substrings. */
  static SAIDPTR
  partition(CSAIDPTR PA,
               SAIDPTR first, SAIDPTR last, SAIDX depth) {
    SAIDPTR a, b;
    SAIDX t;
    for(a = first - 1, b = last;;) {
      for(; (++a < b) && ((PA[*a] + depth) >= (PA[*a + 1] + 1));) { *a = ~*a; }
      for(; (a < --b) && ((PA[*b] + depth) <  (PA[*b + 1] + 1));) { }
      if(b <= a) { break; }
      t = ~*b;
      *b = *a;
      *a = t;
    }
    if(first < a) { *first = ~*first; }
    return a;
  }

  /* Multikey introsort for medium size groups. */
  static void
  mintrosort(CHARPTR T, CSAIDPTR PA,
                SAIDPTR first, SAIDPTR last,
                SAIDX depth) {
#define STACK_SIZE SS_MISORT_STACKSIZE
    std::tuple<SAIDPTR, SAIDPTR, SAIDX, saint_t> stack[STACK_SIZE];
    CHARPTR Td;
    SAIDPTR a, b, c, d, e, f;
    SAIDX s, t;
    saint_t ssize;
    saint_t limit;
    saint_t v, x = 0;

    for(ssize = 0, limit = ilg(last - first);;) {

      if((last - first) <= SS_INSERTIONSORT_THRESHOLD) {
#if 1 < SS_INSERTIONSORT_THRESHOLD
        if(1 < (last - first)) { insertionsort(T, PA, first, last, depth); }
#endif
        if(ssize == 0) return;
        std::tie(first, last, depth, limit) = stack[--ssize];
        continue;
      }

      Td = T + depth;
      if(limit-- == 0) { heapsort(Td, PA, first, last - first); }
      if(limit < 0) {
        for(a = first + 1, v = Td[PA[*first]]; a < last; ++a) {
          if((x = Td[PA[*a]]) != v) {
            if(1 < (a - first)) { break; }
            v = x;
            first = a;
          }
        }
        if(Td[PA[*first] - 1] < v) {
          first = partition(PA, first, a, depth);
        }
        const auto df = a - first, dl = last - a;
        if(df <= dl) {
          if(1 < df) {
            stack[ssize++] = std::make_tuple(a, last, depth, -1);
            last = a, depth += 1, limit = ilg(a - first);
          } else {
            first = a, limit = -1;
          }
        } else {
          if(1 < dl) {
            stack[ssize++] = std::make_tuple(first, a, depth + 1, ilg(a - first));
            first = a, limit = -1;
          } else {
            last = a, depth += 1, limit = ilg(a - first);
          }
        }
        continue;
      }

      /* choose pivot */
      a = pivot(Td, PA, first, last);
      v = Td[PA[*a]];
      swap(*first, *a);

      /* partition */
      for(b = first; (++b < last) && ((x = Td[PA[*b]]) == v);) { }
      if(((a = b) < last) && (x < v)) {
        for(; (++b < last) && ((x = Td[PA[*b]]) <= v);) {
          if(x == v) { swap(*b, *a); ++a; }
        }
      }
      for(c = last; (b < --c) && ((x = Td[PA[*c]]) == v);) { }
      if((b < (d = c)) && (x > v)) {
        for(; (b < --c) && ((x = Td[PA[*c]]) >= v);) {
          if(x == v) { swap(*c, *d); --d; }
        }
      }
      for(; b < c;) {
        swap(*b, *c);
        for(; (++b < c) && ((x = Td[PA[*b]]) <= v);) {
          if(x == v) { swap(*b, *a); ++a; }
        }
        for(; (b < --c) && ((x = Td[PA[*c]]) >= v);) {
          if(x == v) { swap(*c, *d); --d; }
        }
      }

      if(a <= d) {
        c = b - 1;

        if((s = a - first) > (t = b - a)) { s = t; }
        for(e = first, f = b - s; 0 < s; --s, ++e, ++f) { swap(*e, *f); }
        if((s = d - c) > (t = last - d - 1)) { s = t; }
        for(e = b, f = last - s; 0 < s; --s, ++e, ++f) { swap(*e, *f); }

        a = first + (b - a), c = last - (d - c);
        b = (v <= Td[PA[*a] - 1]) ? a : partition(PA, a, c, depth);

        if((a - first) <= (last - c)) {
          if((last - c) <= (c - b)) {
            stack[ssize++] = std::make_tuple(b, c, depth + 1, ilg(c - b));
            stack[ssize++] = std::make_tuple(c, last, depth, limit);
            last = a;
          } else if((a - first) <= (c - b)) {
            stack[ssize++] = std::make_tuple(c, last, depth, limit);
            stack[ssize++] = std::make_tuple(b, c, depth + 1, ilg(c - b));
            last = a;
          } else {
            stack[ssize++] = std::make_tuple(c, last, depth, limit);
            stack[ssize++] = std::make_tuple(first, a, depth, limit);
            first = b, last = c, depth += 1, limit = ilg(c - b);
          }
        } else {
          if((a - first) <= (c - b)) {
            stack[ssize++] = std::make_tuple(b, c, depth + 1, ilg(c - b));
            stack[ssize++] = std::make_tuple(first, a, depth, limit);
            first = c;
          } else if((last - c) <= (c - b)) {
            stack[ssize++] = std::make_tuple(first, a, depth, limit);
            stack[ssize++] = std::make_tuple(b, c, depth + 1, ilg(c - b));
            first = c;
          } else {
            stack[ssize++] = std::make_tuple(first, a, depth, limit);
            stack[ssize++] = std::make_tuple(c, last, depth, limit);
            first = b, last = c, depth += 1, limit = ilg(c - b);
          }
        }
      } else {
        limit += 1;
        if(Td[PA[*first] - 1] < v) {
          first = partition(PA, first, last, depth);
          limit = ilg(last - first);
        }
        depth += 1;
      }
    }
#undef STACK_SIZE
  }

#endif /* (SS_BLOCKSIZE == 0) || (SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE) */


  /*---------------------------------------------------------------------------*/

#if SS_BLOCKSIZE != 0
  template<typename SAIDPTR1, typename SAIDPTR2>
  static void
  blockswap(SAIDPTR1 a, SAIDPTR2 b, SAIDX n) {
    SAIDX t;
    for(; 0 < n; --n, ++a, ++b) {
      t = *a, *a = *b, *b = t;
    }
  }

  static void
  rotate(SAIDPTR first, SAIDPTR middle, SAIDPTR last) {
    SAIDPTR a, b;
    SAIDX t;
    SAIDX l, r;
    l = middle - first, r = last - middle;
    for(; (0 < l) && (0 < r);) {
      if(l == r) { blockswap(first, middle, l); break; }
      if(l < r) {
        a = last - 1, b = middle - 1;
        t = *a;
        do {
          *a-- = *b, *b-- = *a;
          if(b < first) {
            *a = t;
            last = a;
            if((r -= l + 1) <= l) { break; }
            a -= 1, b = middle - 1;
            t = *a;
          }
        } while(1);
      } else {
        a = first, b = middle;
        t = *a;
        do {
          *a++ = *b, *b++ = *a;
          if(last <= b) {
            *a = t;
            first = a + 1;
            if((l -= r + 1) <= r) { break; }
            a += 1, b = middle;
            t = *a;
          }
        } while(1);
      }
    }
  }


  /*---------------------------------------------------------------------------*/

  static void
  inplacemerge(CHARPTR T, CSAIDPTR PA,
                  SAIDPTR first, SAIDPTR middle, SAIDPTR last,
                  SAIDX depth) {
    CSAIDPTR p;
    SAIDPTR a, b;
    SAIDX len, half;
    saint_t q, r;
    saint_t x;

    for(;;) {
      auto v = *(last - 1);
      if(v < 0) { x = 1; p = PA + ~v; }
      else      { x = 0; p = PA +  v; }
      for(a = first, len = middle - first, half = len >> 1, r = -1;
          0 < len;
          len = half, half >>= 1) {
        b = a + half;
        v = *b;
        q = compare(T, PA + (0 <= v ? v : ~v), p, depth);
        if(q < 0) {
          a = b + 1;
          half -= (len & 1) ^ 1;
        } else {
          r = q;
        }
      }
      if(a < middle) {
        if(r == 0) { *a = ~*a; }
        rotate(a, middle, last);
        last -= middle - a;
        middle = a;
        if(first == middle) { break; }
      }
      --last;
      if(x != 0) { while(*--last < 0) { } }
      if(middle == last) { break; }
    }
  }


  /*---------------------------------------------------------------------------*/

  /* Merge-forward with internal buffer. */
  static void
  mergeforward(CHARPTR T, CSAIDPTR PA,
                  SAIDPTR first, SAIDPTR middle, SAIDPTR last,
                  SAIDPTR buf, SAIDX depth) {
    SAIDPTR a, b, c, bufend;
    SAIDX t;
    saint_t r;

    bufend = buf + (middle - first) - 1;
    blockswap(buf, first, middle - first);

    for(t = *(a = first), b = buf, c = middle;;) {
      r = compare(T, PA + *b, PA + *c, depth);
      if(r < 0) {
        do {
          *a = *b; ++a;
          if(bufend <= b) { *bufend = t; return; }
          *b = *a; ++b;
        } while(*b < 0);
      } else if(r > 0) {
        do {
          *a = *c; *c = *++a; ++c;
          if(last <= c) {
            while(b < bufend) { *a = *b; *b = *++a; ++b; }
            *a = *b; *b = t;
            return;
          }
        } while(*c < 0);
      } else {
        *c = ~*c;
        do {
          *a = *b; ++a;
          if(bufend <= b) { *bufend = t; return; }
          *b = *a; ++b;
        } while(*b < 0);

        do {
          *a = *c; *c = *++a; ++c;
          if(last <= c) {
            while(b < bufend) { *a = *b; *b = *++a; ++b; }
            *a = *b; *b = t;
            return;
          }
        } while(*c < 0);
      }
    }
  }

  /* Merge-backward with internal buffer. */
  static void
  mergebackward(CHARPTR T, CSAIDPTR PA,
                   SAIDPTR first, SAIDPTR middle, SAIDPTR last,
                   SAIDPTR buf, SAIDX depth) {
    CSAIDPTR p1, p2;
    SAIDPTR a, b, c, bufend;
    SAIDX t;
    saint_t r;
    saint_t x;
    saint_t v;

    bufend = buf + (last - middle) - 1;
    blockswap(buf, middle, last - middle);

    x = 0;
    v = *bufend;
    if(v< 0) { p1 = PA + ~v; x |= 1; }
    else     { p1 = PA +  v; }
    v = *(middle - 1);
    if(v < 0) { p2 = PA + ~v; x |= 2; }
    else      { p2 = PA +  v; }
    for(t = *(a = last - 1), b = bufend, c = middle - 1;;) {
      r = compare(T, p1, p2, depth);
      if(0 < r) {
        if(x & 1) { do { *a = *b; *b = *--a; --b; } while(*b < 0); x ^= 1; }
        *a = *b; --a;
        if(b <= buf) { *buf = t; break; }
        *b = *a; --b;
        v = *b;
        if(v < 0) { p1 = PA + ~v; x |= 1; }
        else      { p1 = PA +  v; }
      } else if(r < 0) {
        if(x & 2) { do { *a = *c; *c = *--a; --c; } while(*c < 0); x ^= 2; }
        *a = *c; *c = *--a; --c;
        if(c < first) {
          while(buf < b) { *a = *b; *b = *--a; --b; }
          *a = *b, *b = t;
          break;
        }
        v = *c;
        if(v < 0) { p2 = PA + ~v; x |= 2; }
        else      { p2 = PA +  v; }
      } else {
        if(x & 1) { do { *a = *b;  *b = *--a; --b; } while(*b < 0); x ^= 1; }
        *a = ~*b; --a;
        if(b <= buf) { *buf = t; break; }
        *b = *a; --b;
        if(x & 2) { do { *a = *c; *c = *--a; --c; } while(*c < 0); x ^= 2; }
        *a = *c; *c = *--a; --c;
        if(c < first) {
          while(buf < b) { *a = *b; *b = *--a; --b; }
          *a = *b, *b = t;
          break;
        }
        v = *b;
        if(v < 0) { p1 = PA + ~v; x |= 1; }
        else      { p1 = PA +  v; }
        v = *c;
        if(v < 0) { p2 = PA + ~v; x |= 2; }
        else      { p2 = PA +  v; }
      }
    }
  }

  /* D&C based merge. */
  static inline SAIDX getidx(SAIDX x) { return 0 <= x ? x : ~x; }
  static void merge_check(CHARPTR T, CSAIDPTR PA, SAIDX depth, SAIDPTR a, SAIDPTR b, saint_t c) {
    if((c & 1) || ((c & 2) && (compare(T, PA + getidx(*(a - 1)), PA + *a, depth) == 0)))
      *a = ~*a;
    if((c & 4) && (compare(T, PA + getidx(*(b - 1)), PA + *b, depth) == 0))
      *b = ~*b;
  }
  static void
  swapmerge(CHARPTR T, CSAIDPTR PA,
            SAIDPTR first, SAIDPTR middle, SAIDPTR last,
            SAIDPTR buf, SAIDX bufsize, SAIDX depth) {
#define STACK_SIZE SS_SMERGE_STACKSIZE
    std::tuple<SAIDPTR, SAIDPTR, SAIDPTR, saint_t> stack[STACK_SIZE];
    SAIDPTR l, r, lm, rm;
    SAIDX m, len, half;
    saint_t ssize;
    saint_t check, next;

    for(check = 0, ssize = 0;;) {
      if((last - middle) <= bufsize) {
        if((first < middle) && (middle < last)) {
          mergebackward(T, PA, first, middle, last, buf, depth);
        }
        merge_check(T, PA, depth, first, last, check);
        if(ssize == 0) return;
        std::tie(first, middle, last, check) = stack[--ssize];
        continue;
      }

      if((middle - first) <= bufsize) {
        if(first < middle) {
          mergeforward(T, PA, first, middle, last, buf, depth);
        }
        merge_check(T, PA, depth, first, last, check);
        if(ssize == 0) return;
        std::tie(first, middle, last, check) = stack[--ssize];
        continue;
      }

      for(m = 0, len = std::min(middle - first, last - middle), half = len >> 1;
          0 < len;
          len = half, half >>= 1) {
        if(compare(T, PA + getidx(*(middle + m + half)),
                   PA + getidx(*(middle - m - half - 1)), depth) < 0) {
          m += half + 1;
          half -= (len & 1) ^ 1;
        }
      }

      if(0 < m) {
        lm = middle - m, rm = middle + m;
        blockswap(lm, middle, m);
        l = r = middle, next = 0;
        if(rm < last) {
          if(*rm < 0) {
            *rm = ~*rm;
            if(first < lm) { for(; *--l < 0;) { } next |= 4; }
            next |= 1;
          } else if(first < lm) {
            for(; *r < 0; ++r) { }
            next |= 2;
          }
        }

        if((l - first) <= (last - r)) {
          stack[ssize++] = std::make_tuple(r, rm, last, (next & 3) | (check & 4));
          middle = lm, last = l, check = (check & 3) | (next & 4);
        } else {
          if((next & 2) && (r == middle)) { next ^= 6; }
          stack[ssize++] = std::make_tuple(first, lm, l, (check & 3) | (next & 4));
          first = r, middle = rm, check = (next & 3) | (check & 4);
        }
      } else {
        if(compare(T, PA + getidx(*(middle - 1)), PA + *middle, depth) == 0) {
          *middle = ~*middle;
        }
        merge_check(T, PA, depth, first, last, check);
        if(ssize == 0) return;
        std::tie(first, middle, last, check) = stack[--ssize];
      }
    }
#undef STACK_SIZE
  }

#endif /* SS_BLOCKSIZE != 0 */


    /*---------------------------------------------------------------------------*/

    /*- Function -*/

    /* Substring sort */
    static void
      sort(CHARPTR T, CSAIDPTR PA,
      SAIDPTR first, SAIDPTR last,
      SAIDPTR buf, SAIDX bufsize,
      SAIDX depth, SAIDX n, saint_t lastsuffix) {
    SAIDPTR a;
#if SS_BLOCKSIZE != 0
    SAIDPTR b, middle, curbuf;
    SAIDX j, k, curbufsize, limit;
#endif
    SAIDX i;

    if(lastsuffix != 0) { ++first; }

#if SS_BLOCKSIZE == 0
    mintrosort(T, PA, first, last, depth);
#else
    if((bufsize < SS_BLOCKSIZE) &&
      (bufsize < (last - first)) &&
      (bufsize < (limit = isqrt(last - first)))) {
      if(SS_BLOCKSIZE < limit) { limit = SS_BLOCKSIZE; }
      buf = middle = last - limit, bufsize = limit;
    } else {
      middle = last, limit = 0;
    }
    for(a = first, i = 0; SS_BLOCKSIZE < (middle - a); a += SS_BLOCKSIZE, ++i) {
#if SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE
      mintrosort(T, PA, a, a + SS_BLOCKSIZE, depth);
#elif 1 < SS_BLOCKSIZE
      insertionsort(T, PA, a, a + SS_BLOCKSIZE, depth);
#endif
      curbufsize = last - (a + SS_BLOCKSIZE);
      curbuf = a + SS_BLOCKSIZE;
      if(curbufsize <= bufsize) { curbufsize = bufsize, curbuf = buf; }
      for(b = a, k = SS_BLOCKSIZE, j = i; j & 1; b -= k, k <<= 1, j >>= 1) {
        swapmerge(T, PA, b - k, b, b + k, curbuf, curbufsize, depth);
      }
    }
#if SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE
    mintrosort(T, PA, a, middle, depth);
#elif 1 < SS_BLOCKSIZE
    insertionsort(T, PA, a, middle, depth);
#endif
    for(k = SS_BLOCKSIZE; i != 0; k <<= 1, i >>= 1) {
      if(i & 1) {
        swapmerge(T, PA, a - k, a, middle, buf, bufsize, depth);
        a -= k;
      }
    }
    if(limit != 0) {
#if SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE
      mintrosort(T, PA, middle, last, depth);
#elif 1 < SS_BLOCKSIZE
      insertionsort(T, PA, middle, last, depth);
#endif
      inplacemerge(T, PA, first, middle, last, depth);
    }
#endif

    if(lastsuffix != 0) {
      /* Insert last type B* suffix. */
      SAIDX PAi[2]; PAi[0] = PA[*(first - 1)], PAi[1] = n - 2;
      for(a = first, i = *(first - 1);
          (a < last) && ((*a < 0) || (0 < compare(T, &(PAi[0]), PA + *a, depth)));
          ++a) {
        *(a - 1) = *a;
      }
      *(a - 1) = i;
    }
    }
}; // struct ss
} // namespace compactsufsort_imp
#endif /* __SSSORT_IMP_H__ */
