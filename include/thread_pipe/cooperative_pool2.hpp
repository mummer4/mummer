#ifndef __THREAD_PIPE_COOPERATIVE_POOL2_H__
#define __THREAD_PIPE_COOPERATIVE_POOL2_H__

#include <cassert>
#include <thread>
#include <atomic>
#include <chrono>

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

// Call the 'produce' method on an object of type D, with 1 or 2
// arguments, based on arity. Check that calling types are correct.
template<int arity, typename D, typename T>
struct call_produce { };
template<typename D, typename T>
struct call_produce<1, D, T> {
  static bool call(D* self, uint32_t i, T& e) { return self->produce(e); }
};
template<typename D, typename T>
struct call_produce<2, D, T> {
  static bool call(D* self, uint32_t i, T& e) { return self->produce(i, e); }
};

template<typename D, typename T>
class cooperative_pool2 {
public:
  typedef circular_buffer<uint32_t> cbT;
  typedef T                         element_type;

private:
  uint32_t       size_;
  element_type*  elts_;
  cbT            cons_prod_;    // FIFO from Consumers to Producers
  cbT            prod_cons_;    // FIFO from Producers to Consumers
  cbT            tokens_;       // FIFO with producer tokens
  const uint32_t max_producers_;
  uint32_t       done_;         // Number of producer that are done

  // RAII token.
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

public:
  cooperative_pool2(uint32_t max_producers = 1, uint32_t size = 4 * std::thread::hardware_concurrency()) :
    size_(size),
    elts_(new element_type[size_]),
    cons_prod_(size_ + 100),
    prod_cons_(size_ + 100),
    tokens_(max_producers + 1),
    max_producers_(max_producers),
    done_(0)
  {
    // Every element is empty and ready to be filled by the producer
    for(size_t i = 0; i < size_; ++i)
      cons_prod_.enqueue_no_check(i);

    // Every producer token is free
    for(uint32_t i = 0; i < max_producers_; ++i)
      tokens_.enqueue_no_check(i);
  }

  ~cooperative_pool2() { delete [] elts_; }

  uint32_t size() const { return size_; }

  element_type* element_begin() { return elts_; }
  element_type* element_end() { return elts_ + size_; }

  // Contains a filled element or is empty. In which case the producer
  // is done and we should stop processing.
  class job {
    cooperative_pool2& cp_;
    uint32_t          i_;       // Index of element
  public:
    job(cooperative_pool2& cp) : cp_(cp), i_(cp_.get_element()) { }
    ~job() { release(); }

    void release() {
      if(!is_empty()) {
        cp_.cons_prod_.enqueue_no_check(i_);
      }
    }
    bool is_empty() const { return i_ == cbT::guard; }
    operator bool() const { return i_ != cbT::guard; }
    void next() {
      release();
      i_ = cp_.get_element();
    }

    element_type& operator*() { return cp_.elts_[i_]; }
    element_type* operator->() { return &cp_.elts_[i_]; }

  private:
    // Disable copy of job
    job(const job& rhs) { }
    job& operator=(const job& rhs) { }
  };
  friend class job;

  /// STL compliant iterator
  class iterator : public std::iterator<std::input_iterator_tag, element_type> {
    job* j_;
  public:
    iterator() : j_(nullptr) { }
    iterator(cooperative_pool2& cp) : j_(new job(cp)) {
      if(j_->is_empty()) {
        delete j_;
        j_ = nullptr;
      }
    }
    iterator(const iterator& rhs) : j_(rhs.j_) { }

    bool operator==(const iterator& rhs) const { return j_ == rhs.j_; }
    bool operator!=(const iterator& rhs) const { return j_ != rhs.j_; }
    element_type& operator*() { return j_->operator*(); }
    element_type* operator->() { return j_->operator->(); }

    iterator& operator++() {
      j_->next();
      if(j_->is_empty()) {
        delete j_;
        j_ = nullptr;
      }
      return *this;
    }

    iterator operator++(int) {
      iterator res(*this);
      ++*this;
      return res;
    }
  };
  iterator begin() { return iterator(*this); }
  const iterator begin() const { return iterator(*this); }
  const iterator end() const { return iterator(); }


private:
  enum PRODUCER_STATUS { PRODUCER_PRODUCED, PRODUCER_DONE, PRODUCER_EXISTS };
  uint32_t get_element() {
    int iteration = 0;

    while(true) {
      // If less than half full -> try to fill up producer to consumer
      // queue. Disregard return value: in any case will
      // attempt to get an element for ourselves
      if(prod_cons_.fill() < prod_cons_.size() / 2)
        become_producer();

      uint32_t i = prod_cons_.dequeue();
      if(i != cbT::guard)
        return i;

      // Try to become producer
      switch(become_producer()) {
      case PRODUCER_PRODUCED:
        iteration = 0; // Produced. Attempt anew to get an element
        break;
      case PRODUCER_DONE:
        return prod_cons_.dequeue();
      case PRODUCER_EXISTS:
        delay(iteration++); // Already a producer. Wait a bit it adds things to queue
        break;
      }
    }
  }

  PRODUCER_STATUS become_producer() {
    typedef utils::function_traits<decltype(&D::produce)> fun_traits;

    if(prod_cons_.is_closed())
      return PRODUCER_DONE;

    // Mark that we have a produce (myself). If not, return. Token
    // will be release automatically at end of method.
    take_token producer_token(tokens_);
    if(!producer_token.has_token())
      return PRODUCER_EXISTS;

    uint32_t i = cbT::guard;
    try {
      while(true) { // Only way out is if produce method is done (returns true or throw an exception)
        i = cons_prod_.dequeue();
        if(i == cbT::guard)
          return PRODUCER_PRODUCED;

        //        if(static_cast<D*>(this)->produce(producer_token.token_, elts_[i])) // produce returns true if done
        if(call_produce<fun_traits::arity, D, T>::call(static_cast<D*>(this), producer_token.token_, elts_[i]))
          break;

        prod_cons_.enqueue_no_check(i);
      }
    } catch(...) { }       // Threw an exception -> same as being done

    // Producing is done for this producer
    cons_prod_.enqueue_no_check(i);
    producer_token.drop();

    uint32_t is_done = ++done_;
    if(is_done < max_producers_)
      return PRODUCER_PRODUCED;

    prod_cons_.close();
    return PRODUCER_DONE;
  }

  // First 16 operations -> no delay. Then exponential back-off up to a second.
  void delay(int iteration) {
    if(iteration < 16)
      return;
    std::this_thread::sleep_for(std::chrono::milliseconds(1 << std::min(iteration - 16, 10)));
  }
};


// Call the 'consume' method on an object of type D, with 1 or 2
// arguments, based on arity. Check that calling types are correct.
template<int arity, typename D, typename T>
struct call_consume { };
template<typename D, typename T>
struct call_consume<1, D, T> {
  static bool call(D* self, uint32_t i, T& e) { return self->consume(e); }
};
template<typename D, typename T>
struct call_consume<2, D, T> {
  static bool call(D* self, uint32_t i, T& e) { return self->consume(i, e); }
};
template<typename D, typename T>
class write_pool {
public:
  typedef circular_buffer<uint32_t> cbT;
  typedef T                         element_type;

private:
  uint32_t       size_;
  element_type*  elts_;
  cbT            cons_prod_;    // FIFO from Consumers to Producers
  cbT            prod_cons_;    // FIFO from Producers to Consumers
  cbT            tokens_;       // FIFO with consumer tokens
  const uint32_t max_consumers_;
  uint32_t       done_;         // Number of consumer that are done

  // RAII token.
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

public:
  write_pool(uint32_t max_consumers = 1, uint32_t size = 4 * std::thread::hardware_concurrency()) :
    size_(size),
    elts_(new element_type[size_]),
    cons_prod_(size_ + 100),
    prod_cons_(size_ + 100),
    tokens_(max_consumers + 1),
    max_consumers_(max_consumers),
    done_(0)
  {
    // Every element is empty and ready to be filled by the producer
    for(size_t i = 0; i < size_; ++i)
      cons_prod_.enqueue_no_check(i);

    // Every consumer token is free
    for(uint32_t i = 0; i < max_consumers_; ++i)
      tokens_.enqueue_no_check(i);
  }

  ~write_pool() { delete [] elts_; }

  uint32_t size() const { return size_; }
  void close() {
    while(true) {
      switch(become_consumer(true)) {
      case CONSUMER_PRODUCED: // Consumed some, try again until no more
        break;
      case CONSUMER_DONE: // If done -> quit
      case CONSUMER_EXISTS: // If consumer exists, let it do the closing -> quit
        return;
      }
    }
  }

  element_type* element_begin() { return elts_; }
  element_type* element_end() { return elts_ + size_; }

  // Contains a filled element or is empty. In which case the consumer
  // is done and we should stop processing.
  class job {
    write_pool& cp_;
    uint32_t          i_;       // Index of element
  public:
    job(write_pool& cp) : cp_(cp), i_(cp_.get_element()) { }
    ~job() { release(); }

    void release() {
      if(!is_empty()) {
        cp_.prod_cons_.enqueue_no_check(i_);
      }
    }
    bool is_empty() const { return i_ == cbT::guard; }
    operator bool() const { return i_ != cbT::guard; }
    void next() {
      release();
      i_ = cp_.get_element();
    }

    element_type& operator*() { return cp_.elts_[i_]; }
    element_type* operator->() { return &cp_.elts_[i_]; }

  private:
    // Disable copy of job
    job(const job& rhs) { }
    job& operator=(const job& rhs) { }
  };
  friend class job;

  /// STL compliant iterator
  class iterator : public std::iterator<std::output_iterator_tag, element_type> {
    job* j_;
  public:
    iterator() : j_(nullptr) { }
    iterator(write_pool& cp) : j_(new job(cp)) {
      if(j_->is_empty()) {
        delete j_;
        j_ = nullptr;
      }
    }
    iterator(const iterator& rhs) : j_(rhs.j_) { }

    bool operator==(const iterator& rhs) const { return j_ == rhs.j_; }
    bool operator!=(const iterator& rhs) const { return j_ != rhs.j_; }
    element_type& operator*() { return j_->operator*(); }
    element_type* operator->() { return j_->operator->(); }

    iterator& operator++() {
      j_->next();
      if(j_->is_empty()) {
        delete j_;
        j_ = nullptr;
      }
      return *this;
    }

    iterator operator++(int) {
      iterator res(*this);
      ++*this;
      return res;
    }
  };
  iterator begin() { return iterator(*this); }
  const iterator begin() const { return iterator(*this); }
  const iterator end() const { return iterator(); }


private:
  enum CONSUMER_STATUS { CONSUMER_PRODUCED, CONSUMER_DONE, CONSUMER_EXISTS };
  uint32_t get_element() {
    int iteration = 0;

    while(true) {
      // If more than half full -> try to empty consumer to producer
      // queue. Disregard return value: in any case will attempt to
      // get an element for ourselves
      if(prod_cons_.fill() > prod_cons_.size() / 2)
        become_consumer();

      uint32_t i = cons_prod_.dequeue();
      if(i != cbT::guard)
        return i;

      // Try to become consumer
      switch(become_consumer()) {
      case CONSUMER_PRODUCED:
        iteration = 0; // Produced. Attempt anew to get an element
        break;
      case CONSUMER_DONE:
        return cons_prod_.dequeue();
      case CONSUMER_EXISTS:
        delay(iteration++); // Already a consumer. Wait a bit it adds things to queue
        break;
      }
    }
  }

  CONSUMER_STATUS become_consumer(bool close = false) {
    typedef utils::function_traits<decltype(&D::consume)> fun_traits;

    if(prod_cons_.is_closed())
      return CONSUMER_DONE;

    // Mark that we have a consume (myself). If not, return. Token
    // will be release automatically at end of method.
    take_token consumer_token(tokens_);
    if(!consumer_token.has_token())
      return CONSUMER_EXISTS;

    uint32_t i = cbT::guard;
    try {
      while(true) { // Only way out is if consume method is done (returns true or throw an exception)
        i = prod_cons_.dequeue();
        if(i == cbT::guard)
          return !close ? CONSUMER_PRODUCED : CONSUMER_DONE;

        if(call_consume<fun_traits::arity, D, T>::call(static_cast<D*>(this), consumer_token.token_, elts_[i]))
          break;

        cons_prod_.enqueue_no_check(i);
      }
    } catch(...) { }       // Threw an exception -> same as being done

    // Consumming is done for this consumer
    prod_cons_.enqueue_no_check(i);
    consumer_token.drop();

    uint32_t is_done = ++done_;
    if(is_done < max_consumers_)
      return CONSUMER_PRODUCED;

    cons_prod_.close();
    return CONSUMER_DONE;
  }

  // First 16 operations -> no delay. Then exponential back-off up to a second.
  void delay(int iteration) {
    if(iteration < 16)
      return;
    std::this_thread::sleep_for(std::chrono::milliseconds(1 << std::min(iteration - 16, 10)));
  }
};

} // namespace imp
} // namespace thread_pipe

#endif /* __THREAD_PIPE_COOPERATIVE_POOL2_H__ */
