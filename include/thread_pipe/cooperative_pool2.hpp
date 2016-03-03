#ifndef __THREAD_PIPE_COOPERATIVE_POOL2_H__
#define __THREAD_PIPE_COOPERATIVE_POOL2_H__

#include <cassert>
#include <thread>
#include <atomic>
#include <chrono>
#include <vector>
#include <type_traits>

#include <thread_pipe/circular_buffer.hpp>
#include <thread_pipe/traits.hpp>

/// Cooperative pool. Provide a link between many producers and many
/// consumers. It is cooperative in the sense that there is no
/// dedicated threads as producer. When the number of elements in the
/// queue from the producer to the consumer is less than half, then
/// the thread requesting an element attempts to become an additional
/// producer. It stays a producer until the producer to consumer queue
/// is full.
///
/// This class must be subclassed using CRTP. `T` is the type of the
/// element passed around in the queues. The derived class must
/// implement the method `bool produce(uint32_t i, T& e)`. It is
/// called when a thread has become a producer. It must set in `e` the
/// new element, unless there is nothing more to produce. It returns
/// `true` if there is nothing more to produce (and `e` is not used),
/// `false` otherwise.
///
/// The maximum number of producers is specified to the constructor of
/// the class (`max_producers`). The parameter `i` passed to `produce`
/// is in [0, max_producers) and it is guaranteed that at any given
/// time, no two producers have the same `i`.
///
/// The following example will produce the integers `[0, 1000 * max)`,
/// with max producers.
///
/// ~~~{.cc}
/// class sequence : public cooperative_pool<sequence, int> {
///   const uint32_t max_;
///   std::vector    cur_;
///   uint32_t       done_;
/// public:
///   sequence(uint32_t max) : max_(max), cur_(max, 0), done_(0) { }
///   bool produce(uint32_t i, int& e) {
///     int& cur = cur_[i];
///     if(cur < max_) {
///       e = i * max_ + cur++;
///       return false;
///     }
///     return true;
///   }
/// };
/// ~~~
///
/// To access the elements (or the jobs) of the sequence, instantiate
/// a `sequence::job` object and check that it is not empty. If empty,
/// the sequence is over.
///
/// ~~~{.cc}
/// sequence seq; // Sequence, instantiated in main thread
/// // In each consumer thread:
/// while(true) {
///   sequence::job j(seq);
///   if(j.is_empty())
///     break;
///   // Do computation using *j and j->
/// }
/// ~~~

namespace thread_pipe {
namespace imp {

// RAII token.
template<typename cbT>
struct take_token {
  cbT&     tokens_;
  uint32_t token_;
  bool     drop_;

  take_token(cbT& tokens) : tokens_(tokens), token_(tokens.dequeue()), drop_(false) { }
  ~take_token() {
    if(has_token() && !drop_) {
      tokens_.enqueue_no_check(token_);
      //        assert(tokens_.enqueue(token_));
    }
  }
  bool has_token() const { return token_ != cbT::guard; }
  void drop() { drop_ = true; }
};

/// Iterator, almost STL compliant. The iterator MUST be run until
/// reaching the end to not "loose" elements in the
/// pool. Alternatively, one can call release() reach the end.
template<typename Pool, typename Category>
class pool_iterator : public std::iterator<Category, typename Pool::element_type> {
  typedef std::iterator<Category, typename Pool::element_type> super;
  typedef typename Pool::size_type                             size_type;

  size_type m_i;
  Pool*     m_cp;
public:
  typedef typename super::value_type        value_type;
  typedef typename super::difference_type   difference_type;
  typedef typename super::pointer           pointer;
  typedef typename super::reference         reference;
  typedef typename super::iterator_category iterator_category;

  pool_iterator() : m_i(Pool::cbT::guard), m_cp(nullptr) { }
  pool_iterator(Pool& cp) : m_i(cp.get_element()), m_cp(&cp) { }
  pool_iterator(const pool_iterator& rhs) : m_i(rhs.m_i), m_cp(rhs.m_cp) { }
  //  ~pool_iterator() { m_cp->release_element(m_i); }

  bool operator==(const pool_iterator& rhs) const { return m_i == rhs.m_i; }
  bool operator!=(const pool_iterator& rhs) const { return m_i != rhs.m_i; }
  value_type& operator*() { return m_cp->elts_[m_i]; }
  value_type* operator->() { return &m_cp->elts_[m_i]; }

  void release() { m_cp->release_element(m_i); }
  pool_iterator& operator++() {
    m_cp->release_element(m_i);
    m_i = m_cp->get_element();
    return *this;
  }

  pool_iterator operator++(int) {
    pool_iterator res(*this);
    ++*this;
    return res;
  }
};


// Call the 'operator()' method on an object of type D, with 1 or 2
// arguments, based on arity. Check that calling types are correct.
template<int arity, typename D, typename T>
struct delegate { };
template<typename D, typename T>
struct delegate<1, D, T> {
  static bool call(D* self, uint32_t i, T& e) { return (*self)(e); }
};
template<typename D, typename T>
struct delegate<2, D, T> {
  static bool call(D* self, uint32_t i, T& e) { return (*self)(i, e); }
};

// Status of the worker
enum ACTIVE_STATUS { WORKED, DONE, EXISTS };

template<typename T, typename S>
class pool_base {
protected:
  typedef S                          size_type;
  typedef circular_buffer<size_type> cbT;
  typedef T                          element_type;

  size_type      size_;
  element_type*  elts_;
  cbT            cons_prod_;    // FIFO from Consumers to Producers
  cbT            prod_cons_;    // FIFO from Producers to Consumers
  cbT            tokens_;       // FIFO with producer tokens
  const uint32_t max_active_;
  uint32_t       done_;         // Number of producer that are done

  // First 16 operations -> no delay. Then exponential back-off up to a second.
  static void delay(int iteration) {
    if(iteration < 16)
      return;
    std::this_thread::sleep_for(std::chrono::milliseconds(1 << std::min(iteration - 16, 10)));
  }

public:
  pool_base(uint32_t max_active, size_type size)
    : size_(size)
    , elts_(new element_type[size_])
    , cons_prod_(size_ + 100)
    , prod_cons_(size_ + 100)
    , tokens_(max_active + 1)
    , max_active_(max_active)
    , done_(0)
  {
    // Every element is empty and ready to be filled by the producer
    for(size_t i = 0; i < size_; ++i)
      cons_prod_.enqueue_no_check(i);

    // Every producer token is free
    for(uint32_t i = 0; i < max_active_; ++i)
      tokens_.enqueue_no_check(i);
  }

  ~pool_base() { delete [] elts_; }

  size_type size() const { return size_; }

  element_type* element_begin() { return elts_; }
  element_type* element_end() { return elts_ + size_; }
};


template<typename D, typename T>
class producer_pool : public pool_base<T, uint32_t> {
  typedef pool_base<T, uint32_t> super;
public:
  typedef typename super::size_type                             size_type;
  typedef circular_buffer<size_type>                            cbT;
  typedef T                                                     element_type;
  typedef pool_iterator<producer_pool, std::input_iterator_tag> iterator;

  friend iterator;

public:
  producer_pool(uint32_t max_active = 1, size_type size = 4 * std::thread::hardware_concurrency())
    : super(max_active, size)
  { }


  iterator begin() { return iterator(*this); }
  const iterator begin() const { return iterator(*this); }
  const iterator end() const { return iterator(); }


private:
  size_type get_element() {
    int iteration = 0;

    while(true) {
      // If less than half full -> try to fill up producer to consumer
      // queue. Disregard return value: in any case will
      // attempt to get an element for ourselves
      if(super::prod_cons_.fill() < super::size() / 2)
        become_producer();

      size_type i = super::prod_cons_.dequeue();
      if(i != cbT::guard)
        return i;

      // Try to become producer
      switch(become_producer()) {
      case WORKED:
        iteration = 0; // Produced. Attempt anew to get an element
        break;
      case DONE:
        return super::prod_cons_.dequeue();
      case EXISTS:
        super::delay(iteration++); // Already a producer. Wait a bit it adds things to queue
        break;
      }
    }
  }

  void release_element(uint32_t i) {
    if(i != cbT::guard)
      super::cons_prod_.enqueue_no_check(i);
  }

  ACTIVE_STATUS become_producer() {
    typedef utils::function_traits<decltype(&D::operator())> fun_traits;
    typedef delegate<fun_traits::arity, D, T> delegate;

    if(super::prod_cons_.is_closed())
      return DONE;

    // Mark that we have a produce (myself). If not, return. Token
    // will be release automatically at end of method.
    take_token<cbT> producer_token(super::tokens_);
    if(!producer_token.has_token())
      return EXISTS;

    size_type i = cbT::guard;
    try {
      while(true) { // Only way out is if produce method is done (returns true or throw an exception)
        i = super::cons_prod_.dequeue();
        if(i == cbT::guard)
          return WORKED;

        //        if(static_cast<D*>(this)->produce(producer_token.token_, elts_[i])) // produce returns true if done
        if(delegate::call(static_cast<D*>(this), producer_token.token_, super::elts_[i]))
          break;

        super::prod_cons_.enqueue_no_check(i);
      }
    } catch(...) { }       // Threw an exception -> same as being done

    // Producing is done for this producer
    super::    cons_prod_.enqueue_no_check(i);
    producer_token.drop();

    uint32_t is_done = ++super::done_;
    if(is_done < super::max_active_)
      return WORKED;

    super::prod_cons_.close();
    return DONE;
  }
};

template<typename D, typename T>
class consumer_pool : public pool_base<T, uint32_t> {
  typedef pool_base<T, uint32_t> super;
public:
  typedef typename super::size_type                              size_type;
  typedef circular_buffer<size_type>                             cbT;
  typedef T                                                      element_type;
  typedef pool_iterator<consumer_pool, std::output_iterator_tag> iterator;

  friend iterator;

public:
  consumer_pool(uint32_t max_active = 1, size_type size = 4 * std::thread::hardware_concurrency())
    : super(max_active, size)
  { }

  void close() {
    while(true) {
      switch(become_consumer(true)) {
      case WORKED: // Consumed some, try again until no more
        break;
      case DONE: // If done -> quit
      case EXISTS: // If consumer exists, let it do the closing -> quit
        return;
      }
    }
  }

  iterator begin() { return iterator(*this); }
  const iterator begin() const { return iterator(*this); }
  const iterator end() const { return iterator(); }


private:
  size_type get_element() {
    int iteration = 0;

    while(true) {
      // If more than half full -> try to empty consumer to producer
      // queue. Disregard return value: in any case will attempt to
      // get an element for ourselves
      if(super::prod_cons_.fill() > super::size() / 2)
        become_consumer();

      size_type i = super::cons_prod_.dequeue();
      if(i != cbT::guard)
        return i;

      // Try to become consumer
      switch(become_consumer()) {
      case WORKED:
        iteration = 0; // Consumed. Attempt anew to get an element
        break;
      case DONE:
        return super::cons_prod_.dequeue();
      case EXISTS:
        super::delay(iteration++); // Already a consumer. Wait a bit it adds things to queue
        break;
      }
    }
  }

  void release_element(size_type i) {
    if(i != cbT::guard)
      super::prod_cons_.enqueue_no_check(i);
  }

  ACTIVE_STATUS become_consumer(bool close = false) {
    typedef utils::function_traits<decltype(&D::operator())> fun_traits;
    typedef delegate<fun_traits::arity, D, T> delegate;

    if(super::cons_prod_.is_closed())
      return DONE;

    // Mark that we have a consume (myself). If not, return. Token
    // will be release automatically at end of method.
    take_token<cbT> consumer_token(super::tokens_);
    if(!consumer_token.has_token())
      return EXISTS;

    size_type i = cbT::guard;
    try {
      while(true) { // Only way out is if consume method is done (returns true or throw an exception)
        i = super::prod_cons_.dequeue();
        if(i == cbT::guard)
          return !close ? WORKED : DONE;

        if(delegate::call(static_cast<D*>(this), consumer_token.token_, super::elts_[i]))
          break;

        super::cons_prod_.enqueue_no_check(i);
      }
    } catch(...) { }       // Threw an exception -> same as being done

    // Consumming is done for this consumer
    super::prod_cons_.enqueue_no_check(i);
    consumer_token.drop();

    uint32_t is_done = ++super::done_;
    if(is_done < super::max_active_)
      return WORKED;

    super::cons_prod_.close();
    return DONE;
  }
};

template<typename Pipe>
class pipe_input_iterator : public std::iterator<std::input_iterator_tag, typename Pipe::element_type> {
  typedef std::iterator<std::input_iterator_tag, typename Pipe::element_type> super;
  typedef pool_iterator<typename Pipe::pool_type, std::input_iterator_tag>    iterator;

  iterator m_it;
  size_t   m_off;

public:
  typedef typename super::value_type        value_type;
  typedef typename super::difference_type   difference_type;
  typedef typename super::pointer           pointer;
  typedef typename super::reference         reference;
  typedef typename super::iterator_category iterator_category;

  pipe_input_iterator(iterator it) : m_it(it), m_off(0) { }

  bool operator==(const pipe_input_iterator& rhs) const { return m_it == rhs.m_it; }
  bool operator!=(const pipe_input_iterator& rhs) const { return m_it != rhs.m_it; }
  value_type& operator*() { return m_it->elts[m_off]; }
  value_type* operator->() { return &this->operator*(); }

  pipe_input_iterator& operator++() {
    ++m_off;
    if(m_off >= m_it->filled) {
      m_it->filled = 0;
      ++m_it;
      m_off = 0;
    }
    return *this;
  }

  pipe_input_iterator operator++(int) {
    pipe_input_iterator res(*this);
    ++*this;
    return res;
  }
};

template<typename Pipe>
class pipe_output_iterator : public std::iterator<std::output_iterator_tag, typename Pipe::element_type> {
  typedef std::iterator<std::output_iterator_tag, typename Pipe::element_type> super;
  typedef pool_iterator<typename Pipe::pool_type, std::output_iterator_tag>    iterator;

  iterator m_it;

public:
  typedef typename super::value_type        value_type;
  typedef typename super::difference_type   difference_type;
  typedef typename super::pointer           pointer;
  typedef typename super::reference         reference;
  typedef typename super::iterator_category iterator_category;

  pipe_output_iterator(iterator it) : m_it(it) { }

  bool operator==(const pipe_output_iterator& rhs) const { return m_it == rhs.m_it; }
  bool operator!=(const pipe_output_iterator& rhs) const { return m_it != rhs.m_it; }
  value_type& operator*() {
    auto& x = *m_it;
    return x.elts[x.filled];
  }
  value_type* operator->() { return &this->operator*(); }
  void flush() { if(m_it->filled) ++m_it; }
  void done() { ++*this; flush(); m_it.release(); }

  pipe_output_iterator& operator++() {
    ++m_it->filled;
    if(m_it->filled >= m_it->elts.size())
      ++m_it;
    return *this;
  }

  pipe_output_iterator operator++(int) {
    pipe_output_iterator res(*this);
    ++*this;
    return res;
  }
};


template<typename T>
struct group {
  size_t         filled;
  std::vector<T> elts;
};

template<typename D, typename T>
class producer_pipe
  : public producer_pool<producer_pipe<D, T>, group<T>>
{
public:
  typedef producer_pool<producer_pipe, group<T>>  pool_type;
  typedef typename pool_type::size_type           size_type;
  typedef T                                       element_type;
  typedef pipe_input_iterator<producer_pipe>      iterator;

  producer_pipe(uint32_t max_active = 1, uint32_t depth = 10, size_type size = 4 * std::thread::hardware_concurrency())
    : pool_type(max_active, size)
  {
    for(auto it = pool_type::element_begin(); it != pool_type::element_end(); ++it) {
      it->filled = 0;
      it->elts.resize(depth);
    }
  }

  bool operator()(uint32_t i, group<T>& g) {
    typedef utils::function_traits<decltype(&D::operator())> fun_traits;
    typedef delegate<fun_traits::arity, D, T> delegate;

    size_t& f = g.filled;
    for(f = 0; f < g.elts.size(); ++f) {
      if(delegate::call(static_cast<D*>(this), i, g.elts[f]))
        break;
    }
    return f == 0; // Nothing filled -> done
  }

  iterator begin() { return iterator(pool_type::begin()); }
  const iterator begin() const { return iterator(pool_type::begin()); }
  const iterator end() const { return iterator(pool_type::end()); }
};

template<typename D, typename T>
class consumer_pipe
  : public consumer_pool<consumer_pipe<D, T>, group<T>>
{
 public:
  typedef consumer_pool<consumer_pipe, group<T>> pool_type;
  typedef typename pool_type::size_type          size_type;
  typedef T                                      element_type;
  typedef pipe_output_iterator<consumer_pipe>    iterator;

  consumer_pipe(uint32_t max_active = 1, uint32_t depth = 10, size_type size = 4 * std::thread::hardware_concurrency())
    : pool_type(max_active, size)
  {
    for(auto it = pool_type::element_begin(); it != pool_type::element_end(); ++it) {
      it->filled = 0;
      it->elts.resize(depth);
    }
  }

  bool operator()(uint32_t i, group<T>& g) {
    typedef utils::function_traits<decltype(&D::operator())> fun_traits;
    typedef delegate<fun_traits::arity, D, T> delegate;

    for(uint32_t f = 0; f < g.filled; ++f) {
      if(delegate::call(static_cast<D*>(this), i, g.elts[f]))
        return true; // Failed -> done
    }
    g.filled = 0;
    return false;
  }

  iterator begin() { return iterator(pool_type::begin()); }
  const iterator begin() const { return iterator(pool_type::begin()); }
  const iterator end() const { return iterator(pool_type::end()); }
};

} // namespace imp
} // namespace thread_pipe

#endif /* __THREAD_PIPE_COOPERATIVE_POOL2_H__ */
