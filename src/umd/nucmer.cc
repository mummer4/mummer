#include "nucmer.hpp"
#include <essaMEM/sparseSA.hpp>
#include <tigr/mgaps.hh>

namespace mummer {
namespace nucmer {

// class FastaRecordSeq {
//   const char* m_seq;
// public:
//   FastaRecordSeq(const char* seq)
//     : m_seq(seq - 1)
//     , m_len(strlen(seq))
//     , m_id
//   { }
  
// };

void match(const char* reference, const char* query, options_type opts) {
  // std::vector<std::string> descr    = { "ref" };
  // std::vector<long>        startpos = { 0 };
  // auto sa = mummer::sparseSA::create_auto(reference, descr, startpos, opts.min_len, true);
  

  // if(opts.orientation & FORWARD) {
  //   auto append_matches = [&](const mummer::match_t& m) { fwd_matches.push_back({ m.ref + 1, m.query + 1, m.len }); };
  //   switch(opts.match) {
  //   case MUM: sa.findMUM_each(query, opts.min_len, true, append_matches); break;
  //   case MUMREFERENCE: sa.findMAM_each(query, opts.min_len, true, append_matches); break;
  //   case MAXMATCH: sa.findMEM_each(query, opts.min_len, true, append_matches); break;
  //   }
  // }
  // if(opts.orientation & REVERSE) {
  //   auto append_matches = [&](const mummer::match_t& m) { fwd_matches.push_back({ m.ref + 1, m.query + 1, m.len }); };
  //   switch(opts.match) {
  //   case MUM: sa.findMUM_each(query, opts.min_len, false, append_matches); break;
  //   case MUMREFERENCE: sa.findMAM_each(query, opts.min_len, false, append_matches); break;
  //   case MAXMATCH: sa.findMEM_each(query, opts.min_len, false, append_matches); break;
  //   }
  // }

  
}
} // namespace nucmer
} // namespace mummer
