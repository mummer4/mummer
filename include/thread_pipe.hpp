#ifndef __THREAD_PIPE_H__
#define __THREAD_PIPE_H__

#include <iterator>
#include <iostream>
#include <sstream>

#include <thread_pipe/cooperative_pool2.hpp>

namespace thread_pipe {
template<typename D, typename T>
using producer = imp::producer_pipe<D, T>;
template<typename D, typename T>
using consumer = imp::consumer_pipe<D, T>;

// Parse type T (with operator>>) from istream
template<typename T>
class from_istream : public producer<from_istream<T>, T> {
  std::istream& is_;
public:
  from_istream(std::istream& is) : is_(is) { }
  bool operator()(T& e) { return !(is_ >> e); }
};

// Write type T (with operator<<) to ostream
template<typename T>
class to_ostream : public consumer<to_ostream<T>, T> {
  std::ostream& os_;
  const char*   delim_;
public:
  to_ostream(std::ostream& os, const char* delim = "\n") : os_(os), delim_(delim) { }
  bool operator()(T& e) { return !(os_ << e << delim_); }
};

// Iterator behaves like a pointer to an ostream.
struct stringstream_wrapper {
  std::stringstream* p_;
  stringstream_wrapper() : p_(new std::stringstream) { }
  ~stringstream_wrapper() { delete p_; }
  operator std::ostream&() { return *p_; }
  ssize_t tellp() { return p_->tellp(); }
};
template<typename T>
std::ostream& operator<<(stringstream_wrapper& os, const T& x) {
  return *os.p_ << x;
}
class ostream_buffered : public consumer<ostream_buffered, stringstream_wrapper> {
  std::ostream& os_;
public:
  ostream_buffered(std::ostream& os) : os_(os) { }
  bool operator()(stringstream_wrapper& e) {
    bool res = false;
    auto rdbuf = e.p_->rdbuf();
    if(rdbuf->in_avail()) {
      res = !(os_ << rdbuf);
      e.p_->str("");
    }
    return res;
  }
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
  bool operator()(typename super::element_type& e) {
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
