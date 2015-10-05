#ifndef __COMMON_H__
#define __COMMON_H__

#include <istream>
#include <limits>
#include <vector>

inline std::istream& skip_line(std::istream& is) {
  return is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

struct entry {
  int                      size;
  std::vector<std::string> lines;
  void add_line(std::istream& is) {
    if(size >= (int)lines.size())
      lines.push_back("");
    std::getline(is, lines[size]);
    ++size;
  }
};

#endif /* __COMMON_H__ */
