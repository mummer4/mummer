#include "common.hpp"

const char* find_token(uint32_t nb, const char* str) {
  static const char* space = " \t\n";
  for(uint32_t i = 1; i < nb && *str && *str != '\n'; ++i) {
    str += strcspn(str, space);
    str += strspn(str, space);
  }
  return (*str && *str != '\n') ? str : nullptr;
}

size_t append_line(std::istream& is, std::vector<char>& buf, size_t off) {
  if(buf.size() < off + 1024)
    buf.resize(std::max(buf.size() + (buf.size() >> 1), off + 1024));

  is.getline(buf.data() + off, buf.size() - off);
  const size_t r = is.gcount();
  if((is.rdstate() & std::istream::failbit) != 0) {
    is.clear(is.rdstate() ^ std::istream::failbit);
    if(r == 0) return off;
    return append_line(is, buf, off + r);
  }

  if(r > 0)
    buf[off+r-1] = '\n';
  return off + r;
}
