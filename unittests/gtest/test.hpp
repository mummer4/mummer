#ifndef __TEST_H__
#define __TEST_H__

#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include <libgen.h>
#include <string.h>


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

struct prefix_unlink {
  const std::string path;
  bool do_unlink;
  explicit prefix_unlink(const char* s, bool d = true) : path(s), do_unlink(d) { }
  explicit prefix_unlink(const std::string& s, bool d = true) : path(s), do_unlink(d) { }
  ~prefix_unlink() {
    if(do_unlink) {
      std::vector<char> dir_mem(path.c_str(), path.c_str() + path.size() + 1);
      std::vector<char> base_mem(dir_mem);
      const char* dir = dirname(dir_mem.data());
      const char* base = basename(base_mem.data());

      DIR* dirp = opendir(dir);
      EXPECT_NE(nullptr, dirp);
      dirent* ent;
      while((ent = readdir(dirp)) != nullptr) {
        if(strncmp(base, ent->d_name, strlen(base)) == 0)
          unlink(ent->d_name);
      }
      closedir(dirp);
    }
  }
};


// Generate a random sequence of size len
std::string sequence(size_t len);



#endif /* __TEST_H__ */
