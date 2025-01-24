/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __JELLYFISH_CPP_ARRAY_HPP_
#define __JELLYFISH_CPP_ARRAY_HPP_

#include <cstddef>
#include <cstring>
#include <memory>

namespace jellyfish {

/// Fix length array of type T. An element is initialized with the init method.
///   new (this->data() + i) T(
template<typename T>
class cpp_array {
protected:
  // std::pair<T*, ptrdiff_t> data_;
  // std::pair<bool*, ptrdiff_t> init_;
  T* data_;
  bool* init_;
  size_t size_;

public:
  cpp_array(size_t size) :
  data_(static_cast<T*>(::operator new(size * sizeof(T)))),
  init_(new bool[size]),
  size_(size) {
    memset(init_, '\0', sizeof(bool) * size_);
  }

  ~cpp_array() {
    clear();
    ::operator delete(data_); // We don't want to call the destructors: already done by clear()
    delete[] init_;
    // std::return_temporary_buffer(data_.first);
    // std::return_temporary_buffer(init_.first);
  }

  /// Initialize element i with 0 argument
  void init(size_t i) {
    release(i);
    new (data_ + i) T();
    init_[i] = true;
  }

  /// Initialize element i with 1 argument
  template<typename A1>
  void init(size_t i, A1& a1) {
    release(i);
    new (data_ + i) T(a1);
    init_[i] = true;
  }
  template<typename A1>
  void init(size_t i, A1* a1) {
    release(i);
    new (data_ + i) T(a1);
    init_[i] = true;
  }
  /// Initialize element i with 2 arguments
  template<typename A1, typename A2>
  void init(size_t i, A1& a1, A2& a2) {
    release(i);
    new (data_ + i) T(a1, a2);
    init_[i] = true;
  }
  template<typename A1, typename A2>
  void init(size_t i, A1* a1, A2& a2) {
    release(i);
    new (data_ + i) T(a1, a2);
    init_[i] = true;
  }
  template<typename A1, typename A2>
  void init(size_t i, A1& a1, A2* a2) {
    release(i);
    new (data_ + i) T(a1, a2);
    init_[i] = true;
  }
  template<typename A1, typename A2>
  void init(size_t i, A1* a1, A2* a2) {
    release(i);
    new (data_ + i) T(a1, a2);
    init_[i] = true;
  }

  /// Initialize element i with 3 arguments
  template<typename A1, typename A2, typename A3>
  void init(size_t i, A1 a1, A2 a2, A3 a3) {
    release(i);
    new (data_ + i) T(a1, a2, a3);
    init_[i] = true;
  }
  /// Initialize element i with 4 arguments
  template<typename A1, typename A2, typename A3, typename A4>
  void init(size_t i, A1 a1, A2 a2, A3 a3, A4 a4) {
    release(i);
    new (data_ + i) T(a1, a2, a3, a4);
    init_[i] = true;
  }
  /// Initialize element i with 5 arguments
  template<typename A1, typename A2, typename A3, typename A4, typename A5>
  void init(size_t i, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5) {
    release(i);
    new (data_ + i) T(a1, a2, a3, a4, a5);
    init_[i] = true;
  }

  void release(size_t i) {
    if(init_[i]) {
      data_[i].~T();
      init_[i] = false;
    }
  }

  size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  T& operator[](size_t i) { return data_[i]; }
  const T& operator[](size_t i) const { return data_[i]; }
  bool initialized(size_t i) const { return init_[i]; }

  T* begin() { return data_; }
  T* end() { return data_ + size_; }
  const T* begin() const { return data_; }
  const T* end() const { return data_ + size_; }
  const T* cbegin() const { return data_; }
  const T* cend() const { return data_ + size_; }

  T* data() { return data_; }
  const T* data() const { return data_; }

  T& front() { return data_[0]; }
  T& back() { return data_[size_ - 1]; }
  const T& front() const { return data_[0]; }
  const T& back() const { return data_[size_ - 1]; }

  void clear() {
    for(size_t i = 0; i < size_; ++i)
      release(i);
  }
};
} // namespace jellyfish

#endif /* __JELLYFISH_CPP_ARRAY_HPP_ */
