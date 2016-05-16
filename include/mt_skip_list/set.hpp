#ifndef __MULTI_THREAD_SKIP_LIST_SET_HPP__
#define __MULTI_THREAD_SKIP_LIST_SET_HPP__

#include <assert.h>

#include <algorithm>
#include <utility>
#include <functional>
#include <atomic>
#include <fstream>
#include <initializer_list>
#include <stdexcept>

#include <mt_skip_list/common.hpp>

/** Multi-threaded lock free set based on skip list. In multi-threaded
 * mode, it supports only adding an element to the set and testing for
 * membership. Only the member functions of the class
 * multi_thread_skip_list_set which marked "const" are thread
 * safe. Use the method of the class
 * multi_thread_skip_list_set::thread to insert elements. It is
 * constructed from a multi_thread_skip_list_set class.
 */

/* TODO: There is a lot of code duplication with
   skip_list_set.hpp. Can it be reduced while still being readable?
 */

namespace mt_skip_list {

template <typename Key, typename Compare = std::less<Key>, int p_ = 4, typename Random = imp::xor_random>
class set {
protected:
  struct node;

  // Atomic pointer to a node. Default access has relaxed memory order
  struct anode : public std::atomic<node*> {
    operator node*() { return this->load(std::memory_order_relaxed); }
    node* operator=(node* val) {
      this->store(val, std::memory_order_relaxed);
      return val;
    }
  };

  // Node type in the skip list
  struct node {
    Key              k;       // Stored key
    std::atomic<int> height;  // Height of tower. Set to 0 before erasing
    anode            tower[]; // Tower of atomic node* of height <height>
  };
  struct path_node {
    anode* ptr;
    node*  val;
  };

  const int  height_upper_bound_;
  anode*     heads_;
  int        max_height_;
  Compare    comp_;

public:
  typedef Key                                   key_type;
  typedef Key                                   value_type;
  typedef Compare                               key_compare;
  typedef Compare                               value_compare;
  typedef Key&                                  reference;
  typedef const Key&                            const_reference;
  typedef size_t                                size_type;
  typedef ssize_t                               difference_type;
  typedef Key*                                  pointer;
  typedef const Key*                            const_pointer;
  static const int p = p_;

  class node_iterator {
  protected:
    node* item;
    node_iterator(node* item_) : item(item_) { }
    node_iterator(const node_iterator& rhs) : item(rhs.item) { }
    void next() { item = item->tower[0]; std::atomic_thread_fence(std::memory_order_acquire); }
  public:
    bool operator==(const node_iterator& rhs) const { return item == rhs.item; }
    bool operator!=(const node_iterator& rhs) const { return item != rhs.item; }
  };

  class iterator :
    public std::iterator<std::forward_iterator_tag, key_type>,
    public node_iterator {
    friend class set;
    iterator(node* item_) : node_iterator(item_) { }
  public:
    iterator() : node_iterator(nullptr) { }
    iterator(const node_iterator& rhs) : node_iterator(rhs) { }

    iterator& operator=(iterator rhs) {
      std::swap(node_iterator::item, rhs.item);
      return *this;
    }
    reference operator*() { return node_iterator::item->k; }
    pointer operator->() { return &node_iterator::item->k; }
    iterator& operator++() {
      node_iterator::next();
      return *this;
    }
    iterator operator++(int) {
      iterator c(*this);
      node_iterator::next();
      return c;
    }
  };
  class const_iterator :
    public std::iterator<std::forward_iterator_tag, key_type>,
    public node_iterator {
    friend class set;
    const_iterator(node* item_) : node_iterator(item_) { }
  public:
    const_iterator() : node_iterator(nullptr) { }
    const_iterator(const const_iterator& rhs) : node_iterator(rhs) { }
    const_iterator(const iterator& rhs) : node_iterator(rhs) { }

    const_iterator& operator=(node_iterator rhs) {
      swap(node_iterator::item, rhs.item);
      return *this;
    }
    const_reference operator*() { return node_iterator::item->k; }
    const_pointer operator->() { return &node_iterator::item->k; }
    const_iterator& operator++() {
      node_iterator::next();
      return *this;
    }
    const_iterator operator++(int) {
      const_iterator c(*this);
      node_iterator::next();
      return c;
    }
  };


  explicit set(const Compare& comp,
               int max_height,
               int height_upper_bound)
    : height_upper_bound_(height_upper_bound)
    , heads_(new anode[height_upper_bound_])
    , max_height_(max_height)
    , comp_(comp)
  { std::fill_n(heads_, height_upper_bound_, nullptr); }

  explicit set(const Compare& comp = Compare())
    : set(comp, 10, imp::height_bound<p_>::value)
  { }

  template<class InputIterator>
  set(InputIterator first, InputIterator last,
      const Compare& comp = Compare())
    : set(comp, 10, imp::height_bound<p_>::value)
  { insert(first, last); }

  set(const std::initializer_list<value_type>& il,
      const Compare& comp = Compare())
    : set(il.begin(), il.end(), comp)
  { }

  set(const set& rhs)
    : set(rhs.comp_, rhs.max_height_, rhs.height_upper_bound_)
  { insert(rhs.cbegin(), rhs.cend()); }

  set(set&& rhs)
    : height_upper_bound_(rhs.height_upper_bound_)
    , heads_(rhs.heads_)
    , max_height_(rhs.max_height_)
    , comp_(rhs.comp_)
  { rhs.heads_ = nullptr; }

  virtual ~set() {
    if(heads_) {
      clear();
      delete [] heads_;
    }
  }

  set& operator=(const set& rhs) {
    set copy(rhs);
    swap(copy);
    return *this;
  }

  set& operator=(set&& rhs) {
    swap(rhs);
    return *this;
  }

  void swap(set& rhs) {
    if(height_upper_bound_ != rhs.height_upper_bound_)
      throw std::invalid_argument("Incompatible height upper bound");
    std::swap(heads_, rhs.heads_);
    std::swap(max_height_, rhs.max_height_);
    std::swap(comp_, rhs.comp_);
  }

  void clear() {
    node* cnode = heads_[0];
    while(cnode) {
      node* nnode = cnode->tower[0];
      delete cnode;
      cnode = nnode;
    }
    std::fill_n(heads_, height_upper_bound_, nullptr);
  }

  /* The following methods are thread safe.
   */
  size_t size() const {
    size_t res = 0;
    for(node* ptr = heads_[0]; ptr; ptr = ptr->tower[0])
      ++res;
    return res;
  }
  bool empty() const { return heads_[0] == nullptr; }
  size_type max_size() const {
    size_type res = 1;
    size_type pp  = p;
    int       n   = max_height_;

    while(n) {
      if(n & 0x1)
        res *= pp;
      pp *= pp;
      n >>= 1;
    }
    return res;
  }
  iterator begin() { return iterator(heads_[0]); }
  const_iterator begin() const { return const_iterator(heads_[0]); }
  const_iterator cbegin() const { return const_iterator(heads_[0]); }
  iterator end() { return iterator(); }
  const_iterator end() const { return const_iterator(); }
  const_iterator cend() const { return const_iterator(); }

  template<typename T>
  size_type count(const T& x) const { return find_node(x) ? 1 : 0; }
  template<typename T>
  std::pair<iterator, iterator> equal_range(const T& x) const {
    static_assert(std::is_pod<path_node>::value, "Path_node is a POD");

    node* ln = nullptr;
    node* n  = find_node(x, ln);
    return n ? std::make_pair(iterator(n), ++iterator(n)) : std::make_pair(iterator(ln), iterator(ln));
  }

  template<typename T>
  iterator find(const T& x) const {
    return iterator(find_node(x));
  }
  template<typename T>
  iterator lower_bound(const T& x) const {
    static_assert(std::is_pod<path_node>::value, "Path_node is a POD");
    node* ln;
    find_node(x, ln);
    return iterator(ln);
  }
  template<typename T>
  iterator upper_bound(const T& x) const {
    static_assert(std::is_pod<path_node>::value, "Path_node is a POD");
    node* ln;
    node* n = find_node(x, ln);
    return n ? ++iterator(n) : iterator(ln);
  }

  std::pair<iterator, bool> insert(const value_type& x) {
    static_assert(std::is_pod<path_node>::value, "Path_node is a POD");
    path_node path[height_upper_bound_];
    int       aheight;
    node*     n = find_node_path(x, path, aheight);
    if(n) return std::make_pair(iterator(n), false);
    return do_insert(path, new_node(x), aheight);
  }

  std::pair<iterator, bool> insert(value_type&& x) {
    static_assert(std::is_pod<path_node>::value, "Path_node is a POD");
    path_node path[height_upper_bound_];
    int       aheight;
    node*     n = find_node_path(x, path, aheight);
    if(n) return std::make_pair(iterator(n), false);
    return do_insert(path, new_node(std::move(x)), aheight);
  }

  template<class InputIterator>
  void insert(InputIterator first, InputIterator last) {
    for( ; first != last; ++first)
      insert(*first);
  }

  void insert(const std::initializer_list<value_type>& il) {
    insert(il.begin(), il.end());
  }

  template<class... Args>
  std::pair<iterator, bool> emplace(Args&&... args) {
    static_assert(std::is_pod<path_node>::value, "Path_node is a POD");
    path_node path[height_upper_bound_];
    int       aheight;
    node *nnode = emplace_node(std::forward<Args>(args)...);
    node*     n = find_node_path(nnode->k, path, aheight);
    if(n) {
      delete nnode;
      return std::make_pair(iterator(n), false);
    }
    return do_insert(path, nnode, aheight);
  }

  size_type erase(const value_type& x) {
    static_assert(std::is_pod<path_node>::value, "Path_node is a POD");
    path_node path[height_upper_bound_];
    int       aheight;
    node*     n = find_node_path(x, path, aheight);
    if(!n) return 0;

    // Try to disconnect node.
  }

private:
  // Find the path to a node equal to x
  template<typename T>
  node* find_node_path(const T& x, path_node path[], int& aheight) const {
    aheight     = max_height_;
    int   i     = aheight - 1;
    node* nnode = nullptr;

    for( ; i >= 0; --i) {
      nnode = heads_[i];
      std::atomic_thread_fence(std::memory_order_acquire);
      if(nnode && comp_(nnode->k, x))
        break;
      path[i].ptr = &heads_[i];
      path[i].val = nnode;
    }

    node* cnode = nnode;
    for( ; i >= 0; --i) {
      nnode = cnode->tower[i];
      std::atomic_thread_fence(std::memory_order_acquire);
      while(nnode && comp_(nnode->k, x)) {
        cnode = nnode;
        nnode = nnode->tower[i];
        std::atomic_thread_fence(std::memory_order_acquire);
      }
      path[i].ptr = &(cnode->tower[i]);
      path[i].val = nnode;
    }

    return nnode && !comp_(x, nnode->k) ? nnode : nullptr;
  }

  // Find a node equal to x. The path is not recorded.
  template<typename T>
  node* find_node(const T& x) const {
    node* lnode;
    return find_node(x, lnode);
  }

  template<typename T>
  node* find_node(const T& x, node*& nnode) const {
    int i = max_height_ - 1;
    nnode = nullptr;

    for( ; i >= 0; --i) {
      nnode = heads_[i];
      std::atomic_thread_fence(std::memory_order_acquire);
      if(nnode && comp_(nnode->k, x))
         break;
    }

    node* cnode = nnode;
    for( ; i >= 0; --i) {
      nnode = cnode->tower[i];
      std::atomic_thread_fence(std::memory_order_acquire);
      while(nnode && comp_(nnode->k, x)) {
        cnode = nnode;
        nnode = cnode->tower[i];
        std::atomic_thread_fence(std::memory_order_acquire);
      }
    }
    return nnode && !comp_(x, nnode->k) ? nnode : nullptr;
  }

  uint64_t new_seed() {
    // Is that good? Should have our own random generator instead of
    // opening a file and relying on the OS?
    std::ifstream is("/dev/urandom");
    uint64_t res;
    is.read((char*)&res, sizeof(res));
    return res;
  }

  // Allocate a new node. Does raw memory allocation of a node with
  // enough space for the tower. Then in place copy construction of
  // the key from x.
  node* alloc_node() {
    static __thread uint64_t seed = 0;
    if(__builtin_expect(seed == 0, 0))
      seed = new_seed();
    const int height = std::min(height_upper_bound_, imp::random_height<Random, p_>::gen(Random::gen(seed)));
    max_height_ = std::max(max_height_, height);
    node* res   = (node*)operator new(sizeof(node) + height * sizeof(anode));
    res->height.store(height, std::memory_order_relaxed);
    new ((void*)&res->tower) anode[height];
    return res;

  }
  node* new_node(const value_type& x) {
    node* res = alloc_node();
    new ((void*)&res->k) value_type(x);
    return res;
  }

  node* new_node(value_type&& x) {
    node* res = alloc_node();
    new ((void*)&res->k) value_type(std::move(x));
    return res;
  }

  template<class... Args>
  node* emplace_node(Args&&... args) {
    node* res = alloc_node();
    new ((void*)&res->k) value_type(std::forward<Args>(args)...);
    return res;
  }

  std::pair<iterator, bool> do_insert(path_node* path, node* n, const int aheight) {
    for(int i = aheight; i < n->height; ++i) { // in case height of n is higher than height so far
      path[i].ptr = &heads_[i];
      path[i].val = *path[i].ptr;
    }
    for(int i = 0; i < n->height; ++i) {
      node*  oval = path[i].val;
      anode* cptr = path[i].ptr;
      n->tower[i] = oval;
      for(node* cval = oval; !cptr->compare_exchange_weak(cval, n, std::memory_order_release); ) {
        std::atomic_thread_fence(std::memory_order_acquire);
        if(comp_(cval->k, n->k)) {
          cptr = &(cval->tower[i]);
          cval = oval;
        } else if(comp_(n->k, cval->k)) {
          n->tower[i] = cval;
          oval        = cval;
        } else {
          delete n;
          return std::make_pair(iterator(cval), false);
        }
      }
    }
    return std::make_pair(iterator(n), true);
  }
};

template <typename Key, typename Comp, int p_, typename Random>
inline void swap(set<Key, Comp, p_, Random>& x, set<Key, Comp, p_, Random>& y) {
  x.swap(y);
}

} // namespace mt_skip_list

#endif /* __MULTI_THREAD_SKIP_LIST_SET_HPP__ */
