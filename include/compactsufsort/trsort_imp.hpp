#ifndef __TRSORT_IMP_H__
#define __TRSORT_IMP_H__

#include <utility>
#include "divsufsort_private.h"
#include "const_iterator_traits.hpp"

using ::std::swap;

namespace compactsufsort_imp {

template<typename SAIDPTR>
struct tr {
  typedef typename type_traits<SAIDPTR>::SAIDX          SAIDX;
  typedef typename const_iterator_traits<SAIDPTR>::type CSAIDPTR;

  /* Simple insertionsort for small size groups. */
  static void
  insertionsort(CSAIDPTR ISAd, SAIDPTR first, SAIDPTR last) {
    SAIDPTR a, b;
    SAIDX t, r;

    for(a = first + 1; a < last; ++a) {
      for(t = *a, b = a - 1; 0 > (r = ISAd[t] - ISAd[*b]);) {
        do { *(b + 1) = *b; } while((first <= --b) && (*b < 0));
        if(b < first) { break; }
      }
      if(r == 0) { *b = ~*b; }
      *(b + 1) = t;
    }
  }


  /*---------------------------------------------------------------------------*/

  static void
  fixdown(CSAIDPTR ISAd, SAIDPTR SA, SAIDX i, SAIDX size) {
    SAIDX j, k;
    SAIDX v;
    SAIDX c, d, e;

    for(v = SA[i], c = ISAd[v]; (j = 2 * i + 1) < size; SA[i] = SA[k], i = k) {
      d = ISAd[SA[k = j++]];
      if(d < (e = ISAd[SA[j]])) { k = j; d = e; }
      if(d <= c) { break; }
    }
    SA[i] = v;
  }

  /* Simple top-down heapsort. */
  static void
  heapsort(CSAIDPTR ISAd, SAIDPTR SA, SAIDX size) {
    SAIDX i, m;
    SAIDX t;

    m = size;
    if((size % 2) == 0) {
      m--;
      if(ISAd[SA[m / 2]] < ISAd[SA[m]]) { swap(SA[m], SA[m / 2]); }
    }

    for(i = m / 2 - 1; 0 <= i; --i) { fixdown(ISAd, SA, i, m); }
    if((size % 2) == 0) { swap(SA[0], SA[m]); fixdown(ISAd, SA, (SAIDX)0, m); }
    for(i = m - 1; 0 < i; --i) {
      t = SA[0], SA[0] = SA[i];
      fixdown(ISAd, SA, (SAIDX)0, i);
      SA[i] = t;
    }
  }


  /*---------------------------------------------------------------------------*/

  /* Returns the median of three elements. */
  static SAIDPTR
  median3(CSAIDPTR ISAd, SAIDPTR v1, SAIDPTR v2, SAIDPTR v3) {
    if(ISAd[*v1] > ISAd[*v2]) { swap(v1, v2); }
    if(ISAd[*v2] > ISAd[*v3]) {
      if(ISAd[*v1] > ISAd[*v3]) { return v1; }
      else { return v3; }
    }
    return v2;
  }

  /* Returns the median of five elements. */
  static SAIDPTR
  median5(CSAIDPTR ISAd,
          SAIDPTR v1, SAIDPTR v2, SAIDPTR v3, SAIDPTR v4, SAIDPTR v5) {
    if(ISAd[*v2] > ISAd[*v3]) { swap(v2, v3); }
    if(ISAd[*v4] > ISAd[*v5]) { swap(v4, v5); }
    if(ISAd[*v2] > ISAd[*v4]) { swap(v2, v4); swap(v3, v5); }
    if(ISAd[*v1] > ISAd[*v3]) { swap(v1, v3); }
    if(ISAd[*v1] > ISAd[*v4]) { swap(v1, v4); swap(v3, v5); }
    if(ISAd[*v3] > ISAd[*v4]) { return v4; }
    return v3;
  }

  /* Returns the pivot element. */
  static SAIDPTR
  pivot(CSAIDPTR ISAd, SAIDPTR first, SAIDPTR last) {
    SAIDPTR middle;
    SAIDX t;

    t = last - first;
    middle = first + t / 2;

    if(t <= 512) {
      if(t <= 32) {
        return median3(ISAd, first, middle, last - 1);
      } else {
        t >>= 2;
        return median5(ISAd, first, first + t, middle, last - 1 - t, last - 1);
      }
    }
    t >>= 3;
    first  = median3(ISAd, first, first + t, first + (t << 1));
    middle = median3(ISAd, middle - t, middle, middle + t);
    last   = median3(ISAd, last - 1 - (t << 1), last - 1 - t, last - 1);
    return median3(ISAd, first, middle, last);
  }


  /*---------------------------------------------------------------------------*/

  struct trbudget_t {
    SAIDX chance;
    SAIDX remain;
    SAIDX incval;
    SAIDX count;

    void init(SAIDX c, SAIDX i) {
      chance = c;
      remain = incval = i;
    }

    saint_t check(SAIDX size) {
      if(size <= remain) { remain -= size; return 1; }
      if(chance == 0) { count += size; return 0; }
      remain += incval - size;
      chance -= 1;
      return 1;
    }
  };

  /*---------------------------------------------------------------------------*/
  static void
  partition(CSAIDPTR ISAd,
            SAIDPTR first, SAIDPTR middle, SAIDPTR last,
            SAIDPTR *pa, SAIDPTR *pb, SAIDX v) {
    SAIDPTR a, b, c, d, e, f;
    SAIDX t, s;
    SAIDX x = 0;

    for(b = middle - 1; (++b < last) && ((x = ISAd[*b]) == v);) { }
    if(((a = b) < last) && (x < v)) {
      for(; (++b < last) && ((x = ISAd[*b]) <= v);) {
        if(x == v) { swap(*b, *a); ++a; }
      }
    }
    for(c = last; (b < --c) && ((x = ISAd[*c]) == v);) { }
    if((b < (d = c)) && (x > v)) {
      for(; (b < --c) && ((x = ISAd[*c]) >= v);) {
        if(x == v) { swap(*c, *d); --d; }
      }
    }
    for(; b < c;) {
      swap(*b, *c);
      for(; (++b < c) && ((x = ISAd[*b]) <= v);) {
        if(x == v) { swap(*b, *a); ++a; }
      }
      for(; (b < --c) && ((x = ISAd[*c]) >= v);) {
        if(x == v) { swap(*c, *d); --d; }
      }
    }

    if(a <= d) {
      c = b - 1;
      if((s = a - first) > (t = b - a)) { s = t; }
      for(e = first, f = b - s; 0 < s; --s, ++e, ++f) { swap(*e, *f); }
      if((s = d - c) > (t = last - d - 1)) { s = t; }
      for(e = b, f = last - s; 0 < s; --s, ++e, ++f) { swap(*e, *f); }
      first += (b - a), last -= (d - c);
    }
    *pa = first, *pb = last;
  }

  static void
  copy(SAIDPTR ISA, CSAIDPTR SA,
       SAIDPTR first, SAIDPTR a, SAIDPTR b, SAIDPTR last,
       SAIDX depth) {
    /* sort suffixes of middle partition
       by using sorted order of suffixes of left and right partition. */
    SAIDPTR c, d, e;
    SAIDX s, v;

    v = b - SA - 1;
    for(c = first, d = a - 1; c <= d; ++c) {
      if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
        *++d = s;
        ISA[s] = d - SA;
      }
    }
    for(c = last - 1, e = d + 1, d = b; e < d; --c) {
      if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
        *--d = s;
        ISA[s] = d - SA;
      }
    }
  }

  static void
  partialcopy(SAIDPTR ISA, CSAIDPTR SA,
              SAIDPTR first, SAIDPTR a, SAIDPTR b, SAIDPTR last,
              SAIDX depth) {
    SAIDPTR c, d, e;
    SAIDX s, v;
    SAIDX rank, lastrank, newrank = -1;

    v = b - SA - 1;
    lastrank = -1;
    for(c = first, d = a - 1; c <= d; ++c) {
      if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
        *++d = s;
        rank = ISA[s + depth];
        if(lastrank != rank) { lastrank = rank; newrank = d - SA; }
        ISA[s] = newrank;
      }
    }

    lastrank = -1;
    for(e = d; first <= e; --e) {
      rank = ISA[*e];
      if(lastrank != rank) { lastrank = rank; newrank = e - SA; }
      if(newrank != rank) { ISA[*e] = newrank; }
    }

    lastrank = -1;
    for(c = last - 1, e = d + 1, d = b; e < d; --c) {
      if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
        *--d = s;
        rank = ISA[s + depth];
        if(lastrank != rank) { lastrank = rank; newrank = d - SA; }
        ISA[s] = newrank;
      }
    }
  }

  static void
  introsort(SAIDPTR ISA, CSAIDPTR ISAd,
            SAIDPTR SA, SAIDPTR first, SAIDPTR last,
            trbudget_t *budget) {
#define STACK_SIZE TR_STACKSIZE
    std::tuple<CSAIDPTR, SAIDPTR, SAIDPTR, saint_t, saint_t> stack[STACK_SIZE];
    const SAIDX incr = ISAd - ISA;
    saint_t limit; //, next;
    saint_t ssize, trlink = -1;

    for(ssize = 0, limit = ilg(last - first);;) {

      if(limit < 0) {
        if(limit == -1) {
          /* tandem repeat partition */
          SAIDPTR a, b;
          partition(ISAd - incr, first, first, last, &a, &b, last - SA - 1);

          /* update ranks */
          if(a < last) {
            const SAIDX v = a - SA - 1;
            for(SAIDPTR c = first; c < a; ++c) { ISA[*c] = v; }
          }
          if(b < last) {
            const SAIDX v = b - SA - 1;
            for(SAIDPTR c = a; c < b; ++c) { ISA[*c] = v; }
          }

          /* push */
          if(1 < (b - a)) {
            stack[ssize++] = std::make_tuple(CSAIDPTR(), a, b, 0, 0);
            stack[ssize++] = std::make_tuple(ISAd - incr, first, last, -2, trlink);
            trlink = ssize - 2;
          }
          if((a - first) <= (last - b)) {
            if(1 < (a - first)) {
              stack[ssize++] = std::make_tuple(ISAd, b, last, ilg(last - b), trlink);
              last = a, limit = ilg(a - first);
            } else if(1 < (last - b)) {
              first = b, limit = ilg(last - b);
            } else {
              if(ssize == 0) return;
              std::tie(ISAd, first, last, limit, trlink) = stack[--ssize];
            }
          } else {
            if(1 < (last - b)) {
              stack[ssize++] = std::make_tuple(ISAd, first, a, ilg(a - first), trlink);
              first = b, limit = ilg(last - b);
            } else if(1 < (a - first)) {
              last = a, limit = ilg(a - first);
            } else {
              if(ssize == 0) return;
              std::tie(ISAd, first, last, limit, trlink) = stack[--ssize];
            }
          }
        } else if(limit == -2) {
          /* tandem repeat copy */
          SAIDPTR a = std::get<1>(stack[--ssize]), b = std::get<2>(stack[ssize]);
          if(std::get<3>(stack[ssize]) == 0) {
            copy(ISA, SA, first, a, b, last, ISAd - ISA);
          } else {
            if(0 <= trlink) { std::get<3>(stack[trlink]) = -1; }
            partialcopy(ISA, SA, first, a, b, last, ISAd - ISA);
          }
          if(ssize == 0) return;
          std::tie(ISAd, first, last, limit, trlink) = stack[--ssize];
        } else {
          /* sorted partition */
          if(0 <= *first) {
            SAIDPTR a = first;
            do { ISA[*a] = a - SA; } while((++a < last) && (0 <= *a));
            first = a;
          }
          if(first < last) {
            SAIDPTR a = first; do { *a = ~*a; } while(*++a < 0);
            saint_t next = (ISA[*a] != ISAd[*a]) ? ilg(a - first + 1) : -1;
            if(++a < last) {
              const SAIDX v = a - SA - 1;
              for(SAIDPTR b = first; b < a; ++b) { ISA[*b] = v; }
            }

            /* push */
            if(budget->check(a - first)) {
              if((a - first) <= (last - a)) {
                stack[ssize++] = std::make_tuple(ISAd, a, last, -3, trlink);
                ISAd += incr, last = a, limit = next;
              } else {
                if(1 < (last - a)) {
                  stack[ssize++] = std::make_tuple(ISAd + incr, first, a, next, trlink);
                  first = a, limit = -3;
                } else {
                  ISAd += incr, last = a, limit = next;
                }
              }
            } else {
              if(0 <= trlink) { std::get<3>(stack[trlink]) = -1; }
              if(1 < (last - a)) {
                first = a, limit = -3;
              } else {
                if(ssize == 0) return;
                std::tie(ISAd, first, last, limit, trlink) = stack[--ssize];
              }
            }
          } else {
            if(ssize == 0) return;
            std::tie(ISAd, first, last, limit, trlink) = stack[--ssize];
          }
        }
        continue;
      }

      if((last - first) <= TR_INSERTIONSORT_THRESHOLD) {
        insertionsort(ISAd, first, last);
        limit = -3;
        continue;
      }

      if(limit-- == 0) {
        heapsort(ISAd, first, last - first);
        for(SAIDPTR a = last - 1; first < a; ) {
          const SAIDX x = ISAd[*a];
          for(--a; (first <= a) && (ISAd[*a] == x); --a) { *a = ~*a; }
        }
        limit = -3;
        continue;
      }

      /* choose pivot */
      SAIDPTR b, a = pivot(ISAd, first, last);
      swap(*first, *a);
      SAIDX pivot = ISAd[*first];

      /* partition */
      partition(ISAd, first, first + 1, last, &a, &b, pivot);
      if((last - first) != (b - a)) {
        saint_t next = (ISA[*a] != pivot) ? ilg(b - a) : -1;

        /* update ranks */
        const SAIDX v = a - SA - 1;
        for(SAIDPTR c = first; c < a; ++c) { ISA[*c] = v; }
        if(b < last) {
          const SAIDX v = b - SA - 1;
          for(SAIDPTR c = a; c < b; ++c) { ISA[*c] = v; }
        }

        /* push */
        if((1 < (b - a)) && (budget->check(b - a))) {
          if((a - first) <= (last - b)) {
            if((last - b) <= (b - a)) {
              if(1 < (a - first)) {
                stack[ssize++] = std::make_tuple(ISAd + incr, a, b, next, trlink);
                stack[ssize++] = std::make_tuple(ISAd, b, last, limit, trlink);
                last = a;
              } else if(1 < (last - b)) {
                stack[ssize++] = std::make_tuple(ISAd + incr, a, b, next, trlink);
                first = b;
              } else {
                ISAd += incr, first = a, last = b, limit = next;
              }
            } else if((a - first) <= (b - a)) {
              if(1 < (a - first)) {
                stack[ssize++] = std::make_tuple(ISAd, b, last, limit, trlink);
                stack[ssize++] = std::make_tuple(ISAd + incr, a, b, next, trlink);
                last = a;
              } else {
                stack[ssize++] = std::make_tuple(ISAd, b, last, limit, trlink);
                ISAd += incr, first = a, last = b, limit = next;
              }
            } else {
              stack[ssize++] = std::make_tuple(ISAd, b, last, limit, trlink);
              stack[ssize++] = std::make_tuple(ISAd, first, a, limit, trlink);
              ISAd += incr, first = a, last = b, limit = next;
            }
          } else {
            if((a - first) <= (b - a)) {
              if(1 < (last - b)) {
                stack[ssize++] = std::make_tuple(ISAd + incr, a, b, next, trlink);
                stack[ssize++] = std::make_tuple(ISAd, first, a, limit, trlink);
                first = b;
              } else if(1 < (a - first)) {
                stack[ssize++] = std::make_tuple(ISAd + incr, a, b, next, trlink);
                last = a;
              } else {
                ISAd += incr, first = a, last = b, limit = next;
              }
            } else if((last - b) <= (b - a)) {
              if(1 < (last - b)) {
                stack[ssize++] = std::make_tuple(ISAd, first, a, limit, trlink);
                stack[ssize++] = std::make_tuple(ISAd + incr, a, b, next, trlink);
                first = b;
              } else {
                stack[ssize++] = std::make_tuple(ISAd, first, a, limit, trlink);
                ISAd += incr, first = a, last = b, limit = next;
              }
            } else {
              stack[ssize++] = std::make_tuple(ISAd, first, a, limit, trlink);
              stack[ssize++] = std::make_tuple(ISAd, b, last, limit, trlink);
              ISAd += incr, first = a, last = b, limit = next;
            }
          }
        } else {
          if((1 < (b - a)) && (0 <= trlink)) { std::get<3>(stack[trlink]) = -1; }
          if((a - first) <= (last - b)) {
            if(1 < (a - first)) {
              stack[ssize++] = std::make_tuple(ISAd, b, last, limit, trlink);
              last = a;
            } else if(1 < (last - b)) {
              first = b;
            } else {
              if(ssize == 0) return;
              std::tie(ISAd, first, last, limit, trlink) = stack[--ssize];
            }
          } else {
            if(1 < (last - b)) {
              stack[ssize++] = std::make_tuple(ISAd, first, a, limit, trlink);
              first = b;
            } else if(1 < (a - first)) {
              last = a;
            } else {
              if(ssize == 0) return;
              std::tie(ISAd, first, last, limit, trlink) = stack[--ssize];
            }
          }
        }
      } else {
        if(budget->check(last - first)) {
          limit = ilg(last - first), ISAd += incr;
        } else {
          if(0 <= trlink) { std::get<3>(stack[trlink]) = -1; }
          if(ssize == 0) return;
          std::tie(ISAd, first, last, limit, trlink) = stack[--ssize];
        }
      }
    }
#undef STACK_SIZE
  }



  /*---------------------------------------------------------------------------*/

  /*- Function -*/

  /* Tandem repeat sort */
  static void
  sort(SAIDPTR ISA, SAIDPTR SA, SAIDX n, SAIDX depth) {
    SAIDPTR ISAd;
    SAIDPTR first, last;
    trbudget_t budget;
    SAIDX t, skip, unsorted;

    budget.init(ilg(n) * 2 / 3, n);
    for(ISAd = ISA + depth; -n < *SA; ISAd += ISAd - ISA) {
      first = SA;
      skip = 0;
      unsorted = 0;
      do {
        if((t = *first) < 0) { first -= t; skip += t; }
        else {
          if(skip != 0) { *(first + skip) = skip; skip = 0; }
          last = SA + ISA[t] + 1;
          if(1 < (last - first)) {
            budget.count = 0;
            introsort(ISA, ISAd, SA, first, last, &budget);
            if(budget.count != 0) { unsorted += budget.count; }
            else { skip = first - last; }
          } else if((last - first) == 1) {
            skip = -1;
          }
          first = last;
        }
      } while(first < (SA + n));
      if(skip != 0) { *(first + skip) = skip; }
      if(unsorted == 0) { break; }
    }
  }
}; // struct tr
} // namespace compactsufsort_imp

#endif /* __TRSORT_IMP_H__ */
