#ifndef __SAID_TRAITS_H__
#define __SAID_TRAITS_H__

namespace compactsufsort_imp {
template<typename SAIDPTR>
struct type_traits {
  typedef typename std::iterator_traits<SAIDPTR>::value_type value_type;
  typedef typename std::remove_const<value_type>::type       SAIDX;
  typedef typename std::add_const<value_type>::type          CSAIDX;
  //  typedef typename const_iterator_trait<SAIDPTR>::type       CSAIDPTR;
};

template<typename CHARPTR>
struct alphabet_traits {
  typedef typename std::iterator_traits<CHARPTR>::value_type value_type;
  static const size_t                                        size = 1 << (8 * sizeof(value_type));
};

} // namespace compactsufsort_imp

#endif /* __SAID_TRAITS_H__ */
