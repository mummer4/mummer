#ifndef __TEST_H__
#define __TEST_H__

#include <unistd.h>

#include <random>

// Random generator initialized with the correct seed.
extern std::mt19937 rand_gen;


// RAII. Automatically unlink a file when object goes out of scope.
struct file_unlink {
  std::string path;
  bool do_unlink;
  explicit file_unlink(const char* s, bool d = true) : path(s), do_unlink(d) { }
  explicit file_unlink(const std::string& s, bool d = true) : path(s), do_unlink(d) { }
  ~file_unlink() {
    if(do_unlink)
      unlink(path.c_str());
  }
};

// Generate a random sequence of size len
std::string sequence(size_t len);



#endif /* __TEST_H__ */
