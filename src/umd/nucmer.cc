#include <cstring>
#include <fstream>
#include "nucmer.hpp"
#include <essaMEM/sparseSA.hpp>
#include <tigr/mgaps.hh>
#include <tigr/postnuc.hh>


namespace mummer {
namespace nucmer {

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

typedef postnuc::Synteny<FastaRecordSeq> synteny_type;

void match(const char* reference, const char* query, std::vector<postnuc::Alignment>& alignments, options_type opts) {

  std::vector<std::string> descr    = { "ref" };
  std::vector<long>        startpos = { 0 };
  auto sa = mummer::sparseSA::create_auto(reference, descr, startpos, opts.min_len, true);
  sa.construct();
  std::vector<mgaps::Match_t> fwd_matches, bwd_matches;
  FastaRecordSeq Ref(reference), Query(query);
  std::vector<synteny_type> syntenys;
  syntenys.push_back(&Ref);
  synteny_type& synteny = syntenys.front();
  mgaps::ClusterMatches clusterer(opts.fixed_separation, opts.max_separation,
                                  opts.min_output_score, opts.separation_factor,
                                  opts.use_extent);
  mgaps::UnionFind UF;

  if(opts.orientation & FORWARD) {
    auto append_matches = [&](const mummer::match_t& m) { fwd_matches.push_back({ m.ref + 1, m.query + 1, m.len }); };
    switch(opts.match) {
    case MUM: sa.findMUM_each(query, opts.min_len, true, append_matches); break;
    case MUMREFERENCE: sa.findMAM_each(query, opts.min_len, true, append_matches); break;
    case MAXMATCH: sa.findMEM_each(query, opts.min_len, true, append_matches); break;
    }
    clusterer.Cluster_each(fwd_matches.data(), UF, fwd_matches.size(), [&](const mgaps::cluster_type& cluster) {
        postnuc::Cluster cl(postnuc::FORWARD_CHAR);
        for(const auto& m : cluster)
          cl.matches.push_back({ m.Start1, m.Start2, m.Len });
        synteny.clusters.push_back(std::move(cl));
      });
  }

  if(opts.orientation & REVERSE) {
    auto append_matches = [&](const mummer::match_t& m) { bwd_matches.push_back({ m.ref + 1, m.query + 1, m.len }); };
    switch(opts.match) {
    case MUM: sa.findMUM_each(query, opts.min_len, true, append_matches); break;
    case MUMREFERENCE: sa.findMAM_each(query, opts.min_len, true, append_matches); break;
    case MAXMATCH: sa.findMEM_each(query, opts.min_len, true, append_matches); break;
    }
    clusterer.Cluster_each(fwd_matches.data(), UF, bwd_matches.size(), [&](const mgaps::cluster_type& cluster) {
        postnuc::Cluster cl(postnuc::REVERSE_CHAR);
        for(const auto& m : cluster)
          cl.matches.push_back({ m.Start1, m.Start2, m.Len });
        synteny.clusters.push_back(std::move(cl));
      });
  }

  postnuc::merge_syntenys merger(opts.do_delta, opts.do_extend, opts.to_seqend, opts.do_shadows,
                                 opts.break_len, opts.banding, sw_align::NUCLEOTIDE);
  merger.processSyntenys_each(syntenys, Query, [&](std::vector<postnuc::Alignment>&& als, const FastaRecordSeq& Af,
                                                   const FastaRecordSeq& Bf) {
      for(auto& al : als)
        alignments.push_back(std::move(al));
                              });
}
} // namespace nucmer
} // namespace mummer
