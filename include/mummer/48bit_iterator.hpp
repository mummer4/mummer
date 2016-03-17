#ifndef __48BIT_ITERATOR_H__
#define __48BIT_ITERATOR_H__

#include "const_iterator_traits.hpp"

// Define the iterator and its const version
template<typename IDX>
class fortyeight_iterator;
template<typename IDX>
class const_fortyeight_iterator;


namespace fortyeight_iterator_imp {

template<typename IDX>
inline IDX get(const uint32_t* p1, const uint16_t* p2) {
  IDX res = *p1 | ((IDX)*p2 << 32);

  if(std::is_signed<IDX>::value && (res & ((IDX)1 << 47)))
    res |= (((IDX)1 << (sizeof(IDX) * 8 - 48)) - 1) << 48;
  return res;
}

template<typename IDX>
inline void set(uint32_t* p1, uint16_t* p2, const IDX x) {
  *p1 = (uint32_t)x;
  *p2 = (uint64_t)x >> 32;
}

template<typename Derived, typename IDX>
class common {
public:
  typedef typename std::iterator<std::random_access_iterator_tag, IDX>::difference_type difference_type;

  void* raw() const { return static_cast<Derived*>(this)->ptr; }

  Derived& operator=(const Derived& rhs) {
    auto self = static_cast<Derived*>(this);
    self->p1 = rhs.p1;
    self->p2 = rhs.p2;
    return *self;
  }
  Derived& operator=(std::nullptr_t p) {
    auto self = static_cast<Derived*>(this);
    self->p1 = nullptr;
    self->p2 = nullptr;
    return *self;
  }

  IDX operator*() const {
    return get<IDX>(static_cast<const Derived*>(this)->p1,
                    static_cast<const Derived*>(this)->p2);
  }

  bool operator==(const Derived& rhs) const {
    return static_cast<const Derived*>(this)->p1 == rhs.p1;
  }
  bool operator!=(const Derived& rhs) const {
    return static_cast<const Derived*>(this)->p1 != rhs.p1;
  }
  bool operator==(std::nullptr_t) {
    return static_cast<const Derived*>(this)->p1 == nullptr;
  }
  bool operator!=(std::nullptr_t) {
    return static_cast<const Derived*>(this)->p1 != nullptr;
  }
  bool operator<(const Derived& rhs) const {
    return static_cast<const Derived*>(this)->p1 < rhs.p1;
  }
  bool operator>(const Derived& rhs) const {
    return static_cast<const Derived*>(this)->p1 > rhs.p1;
  }
  bool operator>=(const Derived& rhs) const {
    return static_cast<const Derived*>(this)->p1 >= rhs.p1;
  }
  bool operator<=(const Derived& rhs) const {
    return static_cast<const Derived*>(this)->p1 <= rhs.p1;
  }
  Derived& operator++() {
    auto self = static_cast<Derived*>(this);
    ++self->p1;
    ++self->p2;
    return *self;
  }
  Derived operator++(int) {
    Derived res(*static_cast<Derived*>(this));
    ++*this;
    return res;
  }
  Derived& operator--() {
    auto self = static_cast<Derived*>(this);
    --self->p1;
    --self->p2;
    return *self;
  }
  Derived operator--(int) {
    Derived res(*static_cast<Derived*>(this));
    --*this;
    return res;
  }
  Derived& operator+=(difference_type n) {
    auto self = static_cast<Derived*>(this);
    self->p1 += n;
    self->p2 += n;
    return *self;
  }
  Derived operator+(difference_type n) const {
    return Derived(*static_cast<const Derived*>(this)) += n;
  }
  Derived& operator-=(difference_type n) {
    auto self = static_cast<Derived*>(this);
    self->p1 -= n;
    self->p2 -= n;
    return *self;
  }
  Derived operator-(difference_type n) const {
    return Derived(*static_cast<const Derived*>(this)) -= n;
  }
  template<typename DD, typename II>
  difference_type operator-(const common<DD, II>& rhs) const {
    return static_cast<const Derived*>(this)->p1 - static_cast<const DD*>(&rhs)->p1;
  }
  IDX operator[](const difference_type n) const {
    return *(*static_cast<const Derived*>(this) + n);
  }
};

template<typename D, typename I>
bool operator==(std::nullptr_t, const common<D, I>& rhs) {
  return rhs == nullptr;
}
template<typename D, typename I>
bool operator!=(std::nullptr_t, const common<D, I>& rhs) {
  return rhs != nullptr;
}

template<typename D, typename I>
D operator+(typename common<D, I>::difference_type lhs, const common<D, I>& rhs) {
  return rhs + lhs;
}

template<typename IDX>
class setter {
  uint32_t* p1;
  uint16_t* p2;
public:
  typedef fortyeight_iterator<IDX> iterator;
  setter(uint32_t* x, uint16_t* y) : p1(x), p2(y) { }
  operator IDX() const { return get<IDX>(p1, p2); }
  setter& operator=(const IDX x) {
    set<IDX>(p1, p2, x);
    return *this;
  }
  setter& operator=(const setter& rhs) {
    return *this = (IDX)rhs;
  }
  iterator operator&() { return iterator(p1, p2); }
};

template<typename IDX>
void swap(setter<IDX>&& x, setter<IDX>&& y) {
  IDX t = x;
  x = (IDX)y;
  y = t;
}

template<typename Derived, typename IDX>
std::ostream& operator<<(std::ostream& os, const common<Derived, IDX>& p) {
  return os << p.raw();
}

} // namespace fortyeight_iterator_imp

template<typename IDX>
class fortyeight_iterator
  : public std::iterator<std::random_access_iterator_tag, IDX>
  , public fortyeight_iterator_imp::common<fortyeight_iterator<IDX>, IDX>
{
  typedef std::iterator<std::random_access_iterator_tag, IDX> super;
  typedef fortyeight_iterator_imp::setter<IDX>                setter_type;

  friend class const_fortyeight_iterator<IDX>;
  friend class fortyeight_iterator_imp::common<fortyeight_iterator<IDX>, IDX>;
  friend class fortyeight_iterator_imp::common<const_fortyeight_iterator<IDX>, IDX>;

  uint32_t* p1;
  uint16_t* p2;

 public:
  typedef typename super::value_type      value_type;
  typedef typename super::difference_type difference_type;

  fortyeight_iterator() = default;
  fortyeight_iterator(uint32_t* x, uint16_t* y) : p1(x), p2(y) { }
  fortyeight_iterator(const fortyeight_iterator& rhs) : p1(rhs.p1), p2(rhs.p2) { }
  fortyeight_iterator(std::nullptr_t ) : p1(nullptr), p2(nullptr) { }

  setter_type operator*() { return setter_type(p1, p2); }
  setter_type operator[](const difference_type n) const { return *(*this + n); }
};

template<typename IDX>
class const_fortyeight_iterator
  : public std::iterator<std::random_access_iterator_tag, const IDX>
  , public fortyeight_iterator_imp::common<const_fortyeight_iterator<IDX>, IDX>
{
  typedef std::iterator<std::random_access_iterator_tag, const IDX> super;
  const uint32_t* p1;
  const uint16_t* p2;

  friend class fortyeight_iterator<IDX>;
  friend class fortyeight_iterator_imp::common<fortyeight_iterator<IDX>, IDX>;
  friend class fortyeight_iterator_imp::common<const_fortyeight_iterator<IDX>, IDX>;

 public:
  typedef typename super::value_type      value_type;
  typedef typename super::difference_type difference_type;

  const_fortyeight_iterator() = default;
  const_fortyeight_iterator(const uint32_t* x, const uint16_t* y) : p1(x), p2(y) { }
  const_fortyeight_iterator(const const_fortyeight_iterator& rhs) : p1(rhs.p1), p2(rhs.p2) { }
  const_fortyeight_iterator(const fortyeight_iterator<IDX>& rhs) : p1(rhs.p1), p2(rhs.p2) { }
  const_fortyeight_iterator(std::nullptr_t ) : p1(nullptr), p2(nullptr) { }
};

// Type traits
namespace compactsufsort_imp {
template<typename I>
struct const_iterator_traits<fortyeight_iterator<I>> {
  typedef const_fortyeight_iterator<I> type;
};

template<typename I>
struct const_iterator_traits<const_fortyeight_iterator<I>> {
  typedef const_fortyeight_iterator<I> type;
};
} // namespace compactsufsort_imp

#endif /* __48BIT_ITERATOR_H__ */
