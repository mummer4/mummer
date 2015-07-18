#ifndef __NUCMER_H__
#define __NUCMER_H__

#include <vector>
#include <tigr/postnuc.hh>

namespace mummer {
namespace nucmer {
enum match_type { MUM, MUMREFERENCE, MAXMATCH };
enum ori_type { FORWARD = 1, REVERSE = 2, BOTH = 3 };
struct options_type {
  options_type()
    : match(MUMREFERENCE)
    , min_len(20)
    , orientation(BOTH)
    , fixed_separation(5)
    , max_separation(1000)
    , min_output_score(200)
    , separation_factor(0.05)
    , use_extent(false)
    , do_delta(true)
    , do_extend(true)
    , to_seqend(false)
    , do_shadows(false)
  { }

  // Options for mummer
  match_type match;
  int        min_len;
  ori_type   orientation;

  // Options for mgaps
  long   fixed_separation;
  long   max_separation;
  long   min_output_score;
  double separation_factor;
  bool   use_extent;

  // Options for postnuc
  bool do_delta;
  bool do_extend;
  bool to_seqend;
  bool do_shadows;
};

// Match a single DNA sequence <query> against a single DNA sequence <reference>.
void match(const char* reference, const char* query, std::vector<postnuc::Alignment>& alignments,
           options_type opts = options_type());
inline std::vector<postnuc::Alignment> match(const char* reference, const char* query,
                                             options_type opts = options_type()) {
  std::vector<postnuc::Alignment> alignments;
  match(reference, query, alignments, opts);
  return alignments;
}

} // namespace nucmer
} // namespace mummer

#endif /* __NUCMER_H__ */
