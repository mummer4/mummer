#ifndef __THREAD_PIPE_H__
#define __THREAD_PIPE_H__

#include <iterator>
#include <iostream>

#include <thread_pipe/cooperative_pool2.hpp>

namespace thread_pipe {
template<typename D, typename T>
using producer = imp::cooperative_pool2<D, T>;
template<typename D, typename T>
using consumer = imp::write_pool<D, T>;

template<typename T>
class from_istream : public producer<from_istream<T>, T> {
  std::istream& is_;
public:
  from_istream(std::istream& is) : is_(is) { }
  bool produce(T& e) { return !(is_ >> e); }
};

template<typename T>
class to_ostream : public consumer<to_ostream<T>, T> {
  std::ostream& os_;
  const char*   delim_;
public:
  to_ostream(std::ostream& os, const char* delim = "\n") : os_(os), delim_(delim) { }
  bool consume(T& e) { return !(os_ << e << delim_); }
};

template<typename I>
class input_iterator : public producer<input_iterator<I>, typename std::iterator_traits<I>::value_type> {
  I       f_;
  const I l_;
  typedef producer<input_iterator<I>, typename std::iterator_traits<I>::value_type> super;
public:
  input_iterator(I first, I last)
    : f_(first), l_(last)
  { }
  template<typename Container>
  input_iterator(Container& c)
    : input_iterator(c.begin(), c.end())
  { }
  bool produce(typename super::element_type& e) {
    bool end = f_ == l_;
    if(!end) {
      e = *f_;
      ++f_;
    }
    return end;
  }
};
} // namespace thread_pipe

#endif /* __THREAD_PIPE_H__ */
