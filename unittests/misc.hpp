#ifndef __MISC_HPP__
#define __MISC_HPP__

#include <cstring>
#include <stdexcept>
#include <thread>
#include <utility>

template <class Fn, class... Args>
void pdo(unsigned int n, Fn&& fn, Args&&... args) {
  std::vector<std::thread> threads;
  for(unsigned int i = 0; i < n; ++i)
    threads.push_back(std::thread(std::forward<Fn>(fn), std::forward<Args>(args)...));
  for(auto& th : threads)
    th.join();
}

#endif /* __MISC_HPP__ */
