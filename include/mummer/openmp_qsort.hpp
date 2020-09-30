#ifndef __OPENMP_QSORT_H__
#define __OPENMP_QSORT_H__

#include <algorithm>

namespace openmp_qsort_imp {
template<typename Iterator, class Compare>
void openmp_qsort_(Iterator begin, Iterator end, size_t sz, Compare Comp);
} // namespace openmp_qsort_imp

template<typename Iterator, class Compare>
void openmp_qsort(Iterator begin, Iterator end, Compare Comp) {
  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
  static_assert(std::is_same<std::random_access_iterator_tag, iterator_category>::value,
                "openmp_qsort works only with random iterators");

  auto const sz = end - begin;
  if(sz < 1024)
    return std::sort(begin, end, Comp);

//#pragma omp parallel
  {
//#pragma omp single
    {
      openmp_qsort_imp::openmp_qsort_(begin, end, sz, Comp);
    }
  }
}

template<typename Iterator>
void openmp_qsort(Iterator begin, Iterator end) {
  typedef typename std::iterator_traits<Iterator>::value_type Type;
  openmp_qsort(begin, end, std::less<Type>());
}

namespace openmp_qsort_imp {
template<typename Iterator, class Compare>
void openmp_qsort_(Iterator begin, Iterator end, const size_t sz, Compare Comp) {
  typedef typename std::iterator_traits<Iterator>::value_type Type;
  assert((size_t)(end - begin) == sz);
  auto pivot = begin + sz/2;
  auto const pivot_v = *pivot;

  std::swap(*pivot, *(end - 1));
  auto p = std::partition(begin, end, [&](const Type& a) { return Comp(a, pivot_v); });
  std::swap(*p, *(end - 1));

  auto const sz1 = p - begin, sz2 = end - p - 1;
  assert(sz1 >= 0);
  assert(sz2 >= 0);
  assert((size_t)sz2 <= sz);
  assert((size_t)sz1 <= sz);
  assert((size_t)sz1 + (size_t)sz2 + 1 == sz);
  if(sz1 > 1024) {
//#pragma omp task firstprivate(p, sz1)
    openmp_qsort_(begin, p, sz1, Comp);
    if(sz2 > 1024)
      openmp_qsort_(p + 1, end, sz2, Comp);
    else
      std::sort(p + 1, end, Comp);
  } else {
    if(sz2 > 1024)
//#pragma omp task firstprivate(p, sz2)
      openmp_qsort_(p + 1, end, sz2, Comp);
    else
      std::sort(p + 1, end, Comp);
    std::sort(begin, p, Comp);
  }
}

} // namespace openmp_qsort_imp

#endif /* __OPENMP_QSORT_H__ */
