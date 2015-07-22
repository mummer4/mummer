#ifndef __NUCMER_H__
#define __NUCMER_H__

#include <vector>
#include <essaMEM/sparseSA.hpp>
#include <tigr/mgaps.hh>
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
  const options_type            options;

public:
  SequenceAligner(const char* reference, options_type opts = options_type())
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
                            options_type opts = options_type()) {
  SequenceAligner aligner(reference, opts);
  aligner.align(query, alignments);
}

inline std::vector<postnuc::Alignment> align_sequences(const char* reference, const char* query,
                                                       options_type opts = options_type()) {
  SequenceAligner aligner(reference, opts);
  return aligner.align(query);
}

} // namespace nucmer
} // namespace mummer

#endif /* __NUCMER_H__ */
