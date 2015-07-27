#ifndef __NUCMER_H__
#define __NUCMER_H__

#include <vector>
#include <mummer/sparseSA.hpp>
#include <mummer/mgaps.hh>
#include <mummer/postnuc.hh>

namespace mummer {
namespace nucmer {
enum match_type { MUM, MUMREFERENCE, MAXMATCH };
enum ori_type { FORWARD = 1, REVERSE = 2, BOTH = 3 };
struct Options {
  Options()
    : match(MUMREFERENCE)
    , min_len(20)
    , orientation(BOTH)
    , fixed_separation(5)
    , max_separation(90)
    , min_output_score(65)
    , separation_factor(0.12)
    , use_extent(false)
    , do_delta(true)
    , do_extend(true)
    , to_seqend(false)
    , do_shadows(false)
    , break_len(200)
    , banding(0)
  { }

  // Setters corresponding to nucmer.pl switches
  Options& mum() { match = MUM; return *this; }
  Options& mumcand() { return mumreference(); }
  Options& mumreference() { match = MUMREFERENCE; return *this; }
  Options& maxmatch() { match = MAXMATCH; return *this; }
  Options& breaklen(long l) { break_len = l; return *this; }
  Options& banded() { banding = 1; return *this; }
  Options& nobanded() { banding = 0; return *this; }
  Options& mincluster(long m) { min_output_score = m; return *this; }
  Options& diagdiff(long d) { fixed_separation = d; return *this; }
  Options& diagfactor(double f) { separation_factor = f; return *this; }
  Options& extend() { do_extend = true; return *this; }
  Options& noextend() { do_extend = false; return *this; }
  Options& forward() { orientation = FORWARD; return *this; }
  Options& maxgap(long m) { max_separation = m; return *this; }
  Options& minmatch(long m) { min_len = m; return *this; }
  Options& optimize() { to_seqend = false; return *this; }
  Options& nooptimize() { to_seqend = true; return *this; }
  Options& reverse() { orientation = REVERSE; return *this; }
  Options& simplify() { do_shadows = false; return *this; }
  Options& nosimplify() { do_shadows = true; return *this; }

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
  int  break_len;
  int  banding;
};

// FastaRecord information, pointing to an existing string
class FastaRecordSeq {
  const char* m_seq;
  size_t      m_len;
  const std::string m_Id;

public:
  FastaRecordSeq(const char* seq, const char* Id = "")
    : m_seq(seq - 1)
    , m_len(strlen(seq))
    , m_Id(Id)
  { }
  FastaRecordSeq(const char* seq, long int len, const char* Id = "")
    : m_seq(seq - 1)
    , m_len(len)
    , m_Id(Id)
  { }
  const std::string& Id() const { return m_Id; }
  long len() const { return m_len; }
  const char* seq() const { return m_seq; }
};


class SequenceAligner {
  static const std::vector<std::string> descr;
  static const std::vector<long>        startpos;

  const mummer::sparseSA        sa;
  const mgaps::ClusterMatches   clusterer;
  const postnuc::merge_syntenys merger;
  const FastaRecordSeq          Ref;
  const Options                 options;

public:
  SequenceAligner(const char* reference, Options opts = Options())
    : sa(mummer::sparseSA::create_auto(reference, descr, startpos, opts.min_len, true))
    , clusterer(opts.fixed_separation, opts.max_separation,
                opts.min_output_score, opts.separation_factor,
                opts.use_extent)
    , merger(opts.do_delta, opts.do_extend, opts.to_seqend, opts.do_shadows,
             opts.break_len, opts.banding, sw_align::NUCLEOTIDE)
    , Ref(reference)
    , options(opts)
  { }

  // Align a single DNA sequence <query> against a single DNA sequence <reference>.
  void align(const char* query, std::vector<postnuc::Alignment>& alignments);

  inline std::vector<postnuc::Alignment> align(const char* query) {
    std::vector<postnuc::Alignment> alignments;
    align(query, alignments);
    return alignments;
  }
};

inline void align_sequences(const char* reference, const char* query, std::vector<postnuc::Alignment>& alignments,
                            Options opts = Options()) {
  SequenceAligner aligner(reference, opts);
  aligner.align(query, alignments);
}

inline std::vector<postnuc::Alignment> align_sequences(const char* reference, const char* query,
                                                       Options opts = Options()) {
  SequenceAligner aligner(reference, opts);
  return aligner.align(query);
}

} // namespace nucmer
} // namespace mummer

#endif /* __NUCMER_H__ */
