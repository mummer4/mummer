#include <cstring>
#include <fstream>
#include <functional>
#include <unordered_map>
#include <stdexcept>
#include <mummer/nucmer.hpp>
#include <mummer/sparseSA.hpp>
#include <mummer/mgaps.hh>
#include <mummer/postnuc.hh>


namespace mummer {
namespace nucmer {

const std::string FastaRecordSeq::empty_string;
typedef postnuc::Synteny<FastaRecordSeq> synteny_type;

static inline char rc(const char c) {
  switch(c) {
  case 'a': return 't';
  case 'c': return 'g';
  case 'g': return 'c';
  case 't': return 'a';
  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';
  }
  return 'n';
}
void reverse_complement(std::string& s) {
  auto st = s.begin();
  auto en = s.end() - 1;

  for( ; st < en; ++st, --en) {
    const char rc_st = rc(*st);
    *st = rc(*en);
    *en = rc_st;
  }
  if(st == en)
    *st = rc(*st);
}

void SequenceAligner::align(const char* query, size_t query_len, std::vector<postnuc::Alignment>& alignments) {
  std::vector<mgaps::Match_t>        fwd_matches(1), bwd_matches(1);
  FastaRecordSeq                     Query(query, query_len);
  std::vector<synteny_type>          syntenys;
  syntenys.push_back(&Ref);
  synteny_type&               synteny = syntenys.front();
  mgaps::UnionFind                   UF;
  char                               cluster_dir;

  auto append_cluster = [&](const mgaps::cluster_type& cluster) {
    postnuc::Cluster cl(cluster_dir);
    cl.matches.push_back({ cluster[0].Start1, cluster[0].Start2, cluster[0].Len });
    for(size_t i = 1; i < cluster.size(); ++i) {
      const auto& m = cluster[i];
      cl.matches.push_back({ m.Start1 + m.Simple_Adj, m.Start2 + m.Simple_Adj, m.Len - m.Simple_Adj });
    }
    synteny.clusters.push_back(std::move(cl));
  };
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
    case MUM: sa.findMUM_each(rquery.c_str(), query_len, options.min_len, false, append_matches); break;
    case MUMREFERENCE: sa.findMAM_each(rquery.c_str(), query_len, options.min_len, false, append_matches); break;
    case MAXMATCH: sa.findMEM_each(rquery.c_str(), query_len, options.min_len, false, append_matches); break;
    }
    cluster_dir = postnuc::REVERSE_CHAR;
    clusterer.Cluster_each(bwd_matches.data(), UF, bwd_matches.size() - 1, append_cluster);
  }

  merger.processSyntenys_each(syntenys, Query,
                              [&](std::vector<postnuc::Alignment>&& als, const FastaRecordSeq& Af,
                                  const FastaRecordSeq& Bf) { for(auto& al : als) alignments.push_back(std::move(al)); });
}

std::unique_ptr<std::ifstream> sequence_info::open_path(const char* path) {
  std::unique_ptr<std::ifstream> data(new std::ifstream(path));
  if(!data->good())
    throw std::runtime_error(std::string("Unable to open '") + path + "'");
  return data;
}

sequence_info::sequence_info(std::istream& data, size_t chunk_size)  {
  std::string meta, line;

  int c = data.peek();
  if(c != '>')
    throw std::runtime_error(std::string("First character must be a '>', got '") + (char)c + "'");

  while(c != EOF && sequence.size() < chunk_size) {
    std::getline(data, line); // Load one line at a time.

    // Read metadata
    const size_t header_offset = headers.size();
    size_t start = line.find_first_not_of(" \t", 1);
    if(start == std::string::npos)
      start = 0;
    const size_t end = line.find_first_of(" \t", start);
    if(end == std::string::npos)
      headers += line.substr(start);
    else
      headers += line.substr(start, end - start);
    headers += '\0';

    // Read sequence
    sequence += '`';
    const size_t sequence_offset = sequence.size();
    for(c = data.peek(); c != EOF && c != '>'; c = data.peek()) {
      for(c = data.get(); c != EOF && c != '\n'; c = data.get()) { // Copy a line
        if(std::isspace(c)) continue;
        sequence += std::tolower(c);
      // std::getline(data, line);
      // const size_t start = line.find_first_not_of(" ");
      // const size_t end = std::min(line.size(), line.find_last_not_of(" \t"));
      // for(size_t i = start; i <= end; ++i)
      //   sequence += std::tolower(line[i]);
      }
    }

    records.push_back({ sequence_offset, header_offset });
  }
  sequence += '`';
  records.push_back({ sequence.size(), headers.size() });
}

FastaRecordPtr sequence_info::find(size_t pos) const {
  auto rec_it = std::upper_bound(records.cbegin(), records.cend(),
                                 pos, [](size_t pos, const record& b) { return pos < b.seq; });
  return FastaRecordPtr(*this, rec_it - records.cbegin() - 1);
}

} // namespace nucmer
} // namespace mummer
