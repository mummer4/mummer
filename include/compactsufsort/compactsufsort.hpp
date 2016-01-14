#ifndef __COMPACTSUFSORT_H__
#define __COMPACTSUFSORT_H__

#include <cstdint>
#include <utility>
#include <iterator>

typedef uint8_t sauchar_t;
typedef int32_t saint_t;

#include "compactsufsort_imp.hpp"

/*- Prototypes -*/
namespace compactsufsort {
using compactsufsort_imp::type_traits;

/**
 * Constructs the suffix array of a given string.
 * @param T[0..n-1] The input string.
 * @param SA[0..n-1] The output array of suffixes.
 * @param n The length of the given string.
 * @return 0 if no error occurred, -1 or -2 otherwise.
 */
template<typename CHARPTR, typename SAIDPTR, typename SAIDX = typename type_traits<SAIDPTR>::SAIDX>
saint_t
create(CHARPTR T, SAIDPTR SA, SAIDX n) { return compactsufsort_imp::SA<CHARPTR, SAIDPTR>::create(T, SA, n); }

/**
 * Checks the correctness of a given suffix array.
 * @param T[0..n-1] The input string.
 * @param SA[0..n-1] The input suffix array.
 * @param n The length of the given string.
 * @param verbose The verbose mode.
 * @return 0 if no error occurred.
 */
template<typename CHARPTR, typename SAIDPTR, typename SAIDX = typename type_traits<SAIDPTR>::SAIDX>
saint_t
check(CHARPTR T, SAIDPTR SA,
      SAIDX n, saint_t verbose) {
  return compactsufsort_imp::SA<CHARPTR, SAIDPTR>::check(T, SA, n, verbose);
}

/**
 * Search for the pattern P in the string T.
 * @param T[0..Tsize-1] The input string.
 * @param Tsize The length of the given string.
 * @param SA[0..Tsize-1] The input suffix array.
 * @param SAsize The length of the suffix array.
 * @param P[0..Psize-1] The input pattern string.
 * @param Psize The length of the given pattern string.
 * @return The count of matches and output index if no error occurred. -1 otherwise.
 */
template<typename CHARPTR, typename CSAIDPTR, typename SAIDX = typename type_traits<CSAIDPTR>::SAIDX>
std::pair<SAIDX, SAIDX>
search(CHARPTR T, SAIDX Tsize, CSAIDPTR SA, SAIDX SAsize,
       CHARPTR P, SAIDX Psize) {
  return compactsufsort_imp::SA<CHARPTR, CSAIDPTR>::search(T, Tsize, SA, SAsize, P, Psize);
}

/**
 * Search for the character c in the string T.
 * @param T[0..Tsize-1] The input string.
 * @param Tsize The length of the given string.
 * @param c The input character.
 * @param SA[0..Tsize-1] The input suffix array.
 * @return The count of matches if no error occurred, -1 otherwise.
 */
template<typename CHARPTR, typename CSAIDPTR, typename SAIDX = typename type_traits<CSAIDPTR>::SAIDX>
SAIDX
search(CHARPTR T, SAIDX Tsize, CSAIDPTR SA, SAIDX SAsize,
       saint_t c) {
  return compactsufsort_imp::SA<CHARPTR, CSAIDPTR>::search(T, Tsize, SA, SAsize, c);
}

} // namespace compactsufsort


#endif /* __COMPACTSUFSORT_H__ */
