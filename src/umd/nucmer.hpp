#ifndef __NUCMER_H__
#define __NUCMER_H__

namespace mummer {
namespace nucmer {
enum match_type { MUM, MUMREFERENCE, MAXMATCH };
enum ori_type { FORWARD = 1, REVERSE = 2, BOTH = 3 };
struct options_type {
  options_type()
    : match(MUMREFERENCE)
    , min_len(20)
    , orientation(BOTH)
  { }

  match_type match;
  int        min_len;
  ori_type   orientation;
};

void match(const char* reference, const char* query, options_type opts = options_type());

} // namespace nucmer
} // namespace mummer

#endif /* __NUCMER_H__ */
