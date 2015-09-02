#ifndef __NUCMER_H__
#define __NUCMER_H__

#include <vector>
#include <forward_list>
#include <thread>
#include <mummer/sparseSA.hpp>
#include <mummer/mgaps.hh>
#include <mummer/postnuc.hh>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>

namespace mummer {
namespace nucmer {
void reverse_complement(std::string& s);

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

// FastaRecord information, pointing to an existing string. Meant to
// be used with postnuc code, which is 1 based. I.e., the first base
// is in m_seq[1]. m_seq[0] contains a separator (probably, don't
// count on it). The addressable part is [1, m_len].
class FastaRecordSeq {
  static const std::string empty_string;
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

///////////////////////////////////////////
// Align two sequences given as strings. //
///////////////////////////////////////////
class SequenceAligner {
  const mummer::sparseSA         sa;
  const mgaps::ClusterMatches    clusterer;
  const postnuc::merge_syntenys  merger;
  const FastaRecordSeq           Ref;
  const Options                  options;

public:
  SequenceAligner(const char* reference, size_t reference_len, const Options opts = Options())
    : sa(mummer::sparseSA::create_auto(reference, reference_len, opts.min_len, true))
    , clusterer(opts.fixed_separation, opts.max_separation,
                opts.min_output_score, opts.separation_factor,
                opts.use_extent)
    , merger(opts.do_delta, opts.do_extend, opts.to_seqend, opts.do_shadows,
             opts.break_len, opts.banding, sw_align::NUCLEOTIDE)
    , Ref(reference)
    , options(opts)
  { }
  explicit SequenceAligner(const std::string& reference, const Options opts = Options())
    : SequenceAligner(reference.c_str(), reference.length(), opts)
  { }

  // Align a single DNA sequence <query> against a single DNA sequence <reference>.
  void align(const char* query, size_t query_len, std::vector<postnuc::Alignment>& alignments);

  inline std::vector<postnuc::Alignment> align(const char* query, size_t query_len) {
    std::vector<postnuc::Alignment> alignments;
    align(query, query_len, alignments);
    return alignments;
  }
};

inline void align_sequences(const char* reference, size_t reference_len,
                            const char* query,  size_t query_len,
                            std::vector<postnuc::Alignment>& alignments,
                            Options opts = Options()) {
  SequenceAligner aligner(reference, reference_len, opts);
  aligner.align(query, query_len, alignments);
}

inline std::vector<postnuc::Alignment> align_sequences(const char* reference, size_t reference_len,
                                                       const char* query,  size_t query_len,
                                                       const Options opts = Options()) {
  SequenceAligner aligner(reference, reference_len, opts);
  return aligner.align(query, query_len);
}


// Sequence concatenated and headers
class FastaRecordPtr;
struct sequence_info {
  struct record { size_t seq, header; };
  std::vector<record> records;
  std::string         sequence;
  std::string         headers;

  // Load from a file
  explicit sequence_info(const char* path);
  sequence_info(const sequence_info& rhs) = delete;
  sequence_info& operator=(const sequence_info& rhs) = delete;
  // Return the FastaRecordPtr corresponding to the sequence containing position pos
  FastaRecordPtr find(size_t pos) const;
};

class FastaRecordPtr {
  const sequence_info& m_info;
  const size_t         m_id;
public:
  FastaRecordPtr(const sequence_info& info, size_t id)
    : m_info(info), m_id(id)
  { assert(m_id < m_info.records.size()); }

  long len() const {
    assert(m_id < m_info.records.size());
    return m_info.records[m_id + 1].seq - m_info.records[m_id].seq - 1;
  }
  const char* seq() const {
    assert(m_id < m_info.records.size());
    return m_info.sequence.c_str() + m_info.records[m_id].seq - 1;
  }
  size_t seq_offset() const {
    assert(m_id < m_info.records.size());
    return m_info.records[m_id].seq;
  }
  const char* Id() const {
    assert(m_id < m_info.records.size());
    return m_info.headers.c_str() + m_info.records[m_id].header;
  }
  bool operator==(const FastaRecordPtr& rhs) const { return m_id == rhs.m_id; }
};


// //////////////////////////////////////////////////////////////////////////////
// // Align many sequences given as sequence files (in fasta or fastq format). //
// //////////////////////////////////////////////////////////////////////////////
class FileAligner {
  const sequence_info           reference_info;
  const mummer::sparseSA        sa;
  const mgaps::ClusterMatches   clusterer;
  const postnuc::merge_syntenys merger;
  const Options                 options;

public:
  FileAligner(const char* reference_path, Options opts = Options())
    : reference_info(reference_path)
    , sa(mummer::sparseSA::create_auto(reference_info.sequence.c_str(), reference_info.sequence.length(), opts.min_len, true))
    , clusterer(opts.fixed_separation, opts.max_separation,
                opts.min_output_score, opts.separation_factor,
                opts.use_extent)
    , merger(opts.do_delta, opts.do_extend, opts.to_seqend, opts.do_shadows,
             opts.break_len, opts.banding, sw_align::NUCLEOTIDE)
    , options(opts)
  { }

  // Align the sequence query against the references
  template<typename AlignmentOut>
  void align(const char* query, size_t query_len, AlignmentOut alignments) const;
  template<typename AlignmentOut>
  void align(const std::string& query, AlignmentOut alignments) const {
    align(query.c_str(), query.length(), alignments);
  }

  // Aligne the sequences in the file query against the references, in
  // parallel
  template<typename AlignmentOut>
  void align_file(const char* query_path, AlignmentOut alignments, unsigned int threads = std::thread::hardware_concurrency()) const;

  typedef jellyfish::stream_manager<const char**>          stream_manager;
  typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;
  template<typename AlignmentOut>
  static void trampoline_align_file(const FileAligner* self, sequence_parser* parser, AlignmentOut alignments) {
    self->thread_align_file(*parser, alignments);
  }
  template<typename AlignmentOut>
  void thread_align_file(sequence_parser& parser, AlignmentOut alignments) const;
};

//
// Implementation of templated methods
//
template<typename AlignmentOut>
void FileAligner::align(const char* query, size_t query_len, AlignmentOut alignments) const {
  typedef postnuc::Synteny<FastaRecordPtr> synteny_type;
  std::vector<mgaps::Match_t>       fwd_matches(1), bwd_matches(1);
  FastaRecordSeq                    Query(query, query_len);
  std::vector<synteny_type>         syntenys;
  std::forward_list<FastaRecordPtr> records;
  mgaps::UnionFind                  UF;
  char                              cluster_dir;

  auto append_cluster = [&](const mgaps::cluster_type& cluster) {
    postnuc::Cluster cl(cluster_dir);
    // Re-map the reference coordinate back to its original sequence
    const auto record  = reference_info.find(cluster[0].Start1);
    const long offset  = record.seq_offset();
    auto       synteny = syntenys.rbegin();
    auto       it      = records.cbegin();
    for( ; it != records.cend(); ++it, ++synteny)
      if(*it == record)
        break;
    if(it == records.cend()) {
      assert(synteny == syntenys.rend());
      records.push_front(record);
      syntenys.push_back(synteny_type(&records.front()));
      synteny = syntenys.rbegin();
    }

    cl.matches.push_back({ cluster[0].Start1 - offset, cluster[0].Start2, cluster[0].Len });
    for(size_t i = 1; i < cluster.size(); ++i) {
      const auto& m = cluster[i];
      cl.matches.push_back({ m.Start1 + m.Simple_Adj - offset, m.Start2 + m.Simple_Adj, m.Len - m.Simple_Adj });
    }
    synteny->clusters.push_back(std::move(cl));
  };

  assert(fwd_matches.size() == 1);
  assert(bwd_matches.size() == 1);
  if(options.orientation & FORWARD) {
    auto append_matches = [&](const mummer::match_t& m) { fwd_matches.push_back({ m.ref + 1, m.query + 1, m.len }); };
    switch(options.match) {
    case MUM: sa.findMUM_each(query, query_len, options.min_len, false, append_matches); break;
    case MUMREFERENCE: sa.findMAM_each(query, query_len, options.min_len, false, append_matches); break;
    case MAXMATCH: sa.findMEM_each(query, query_len, options.min_len, false, append_matches); break;
    }
    cluster_dir = postnuc::FORWARD_CHAR;
    clusterer.Cluster_each(fwd_matches.data(), UF, fwd_matches.size() - 1, append_cluster);
  }

  if(options.orientation & REVERSE) {
    std::string rquery(query, query_len);
    reverse_complement(rquery);
    auto append_matches = [&](const mummer::match_t& m) { bwd_matches.push_back({ m.ref + 1, m.query + 1, m.len }); };
    switch(options.match) {
    case MUM: sa.findMUM_each(rquery, options.min_len, false, append_matches); break;
    case MUMREFERENCE: sa.findMAM_each(rquery, options.min_len, false, append_matches); break;
    case MAXMATCH: sa.findMEM_each(rquery, options.min_len, false, append_matches); break;
    }
    cluster_dir = postnuc::REVERSE_CHAR;
    clusterer.Cluster_each(bwd_matches.data(), UF, bwd_matches.size() - 1, append_cluster);
  }

  merger.processSyntenys_each(syntenys, Query, alignments);
}

template<typename AlignmentOut>
void FileAligner::align_file(const char* query_path, AlignmentOut alignments, unsigned int nb_threads) const {
  stream_manager  streams(&query_path, &query_path + 1);
  sequence_parser parser(4 * nb_threads, 10, 1, streams);

  std::vector<std::thread> threads;
  for(unsigned int i = 0; i < nb_threads; ++i) {
    threads.push_back(std::thread(trampoline_align_file<AlignmentOut>, this, &parser, alignments));
  }
  for(auto& th : threads)
    th.join();
}

template<typename AlignmentOut>
void FileAligner::thread_align_file(sequence_parser& parser, AlignmentOut alignments) const {
  typedef postnuc::Synteny<FastaRecordPtr> synteny_type;
  std::vector<mgaps::Match_t>       fwd_matches(1), bwd_matches(1);
  std::vector<synteny_type>         syntenys;
  std::forward_list<FastaRecordPtr> records;
  mgaps::UnionFind                  UF;
  char                              cluster_dir;

  auto append_cluster = [&](const mgaps::cluster_type& cluster) {
    postnuc::Cluster cl(cluster_dir);
    // Re-map the reference coordinate back to its original sequence
    const auto record  = reference_info.find(cluster[0].Start1);
    const long offset  = record.seq_offset();
    auto       synteny = syntenys.rbegin();
    auto       it      = records.cbegin();
    for( ; it != records.cend(); ++it, ++synteny)
      if(*it == record)
        break;
    if(it == records.cend()) {
      assert(synteny == syntenys.rend());
      records.push_front(record);
      syntenys.push_back(synteny_type(&records.front()));
      synteny = syntenys.rbegin();
    }

    cl.matches.push_back({ cluster[0].Start1 - offset, cluster[0].Start2, cluster[0].Len });
    for(size_t i = 1; i < cluster.size(); ++i) {
      const auto& m = cluster[i];
      cl.matches.push_back({ m.Start1 + m.Simple_Adj - offset, m.Start2 + m.Simple_Adj, m.Len - m.Simple_Adj });
    }
    synteny->clusters.push_back(std::move(cl));
  };

  while(true) {
    sequence_parser::job j(parser);
    if(j.is_empty()) break;

    for(size_t i = 0; i < j->nb_filled; ++i) {
      for(char& c : j->data[i].seq)
        c = std::tolower(c);
      FastaRecordSeq Query(j->data[i].seq.c_str(), j->data[i].seq.length(), j->data[i].header.c_str());
      assert(fwd_matches.size() == 1);
      assert(bwd_matches.size() == 1);
      assert(syntenys.empty());
      if(options.orientation & FORWARD) {
        auto append_matches = [&](const mummer::match_t& m) { fwd_matches.push_back({ m.ref + 1, m.query + 1, m.len }); };
        switch(options.match) {
        case MUM: sa.findMUM_each(Query.seq() + 1, Query.len(), options.min_len, false, append_matches); break;
        case MUMREFERENCE: sa.findMAM_each(Query.seq() + 1, Query.len(), options.min_len, false, append_matches); break;
        case MAXMATCH: sa.findMEM_each(Query.seq() + 1, Query.len(), options.min_len, false, append_matches); break;
        }
        cluster_dir = postnuc::FORWARD_CHAR;
        clusterer.Cluster_each(fwd_matches.data(), UF, fwd_matches.size() - 1, append_cluster);
      }

      if(options.orientation & REVERSE) {
        std::string rquery(Query.seq() + 1, Query.len());
        reverse_complement(rquery);
        auto append_matches = [&](const mummer::match_t& m) { bwd_matches.push_back({ m.ref + 1, m.query + 1, m.len }); };
        switch(options.match) {
        case MUM: sa.findMUM_each(rquery, options.min_len, false, append_matches); break;
        case MUMREFERENCE: sa.findMAM_each(rquery, options.min_len, false, append_matches); break;
        case MAXMATCH: sa.findMEM_each(rquery, options.min_len, false, append_matches); break;
        }
        cluster_dir = postnuc::REVERSE_CHAR;
        clusterer.Cluster_each(bwd_matches.data(), UF, bwd_matches.size() - 1, append_cluster);
      }
      merger.processSyntenys_each(syntenys, Query, alignments);
    }
  }

}

} // namespace nucmer
} // namespace mummer

#endif /* __NUCMER_H__ */
