#ifndef __MUMMER_POSTNUC_H__
#define __MUMMER_POSTNUC_H__

#include <string>
#include <vector>
#include <iostream>
#include <limits>
#include <cstring>
#include <memory>
#include <iomanip>
#include <atomic>

#include "tigrinc.hh"
#include "sw_align.hh"


namespace mummer {
namespace postnuc {
static const signed char FORWARD_CHAR = 1;
static const signed char REVERSE_CHAR = -1;

//------------------------------------------------------ Type Definitions ----//
enum LineType
//-- The type of input line from <stdin>
  {
    NO_LINE, HEADER_LINE, MATCH_LINE
  };


struct Match
//-- An exact match between two sequences A and B
{
  long int sA, sB, len;      // start coordinate in A, in B and the length
};


struct Cluster
//-- An ordered list of matches between two sequences A and B
{
  bool          wasFused;       // have the cluster matches been fused yet?
  signed char   dirB;           // the query sequence direction
                                //      FORWARD_CHAR or REVERSE_CHAR
  std::vector<Match> matches;        // the ordered set of matches in the cluster
  Cluster() = default;
  Cluster(char dir) : wasFused(false), dirB(dir) { }
};

// A FastaRecord must respond to seq(), len() and Id().
template<typename FastaRecord>
struct Synteny
//-- An ordered list of clusters between two sequences A and B (B is fixed)
{
  const FastaRecord*   AfP;     // a pointer to the reference sequence record
  std::vector<Cluster> clusters; // the ordered set of clusters between A and B
  Synteny() = default;
  Synteny(const FastaRecord* Af) : AfP(Af) { }
  bool operator<(const Synteny& rhs) const { return *AfP < *rhs.AfP; }
};

struct Alignment
//-- An alignment object between two sequences A and B
{
  signed char      dirB;        // the query sequence direction
  long int         sA, sB, eA, eB; // the start in A, B and the end in A, B
  std::vector<long int> delta;       // the delta values, with NO zero at the end
  long int         deltaApos;   // sum of abs(deltas) - #of negative deltas
                                //      trust me, it is a very helpful value
  long int         Errors, SimErrors, NonAlphas; // errors, similarity errors, nonalphas

  Alignment(const Match& m, const char dir)
    : dirB(dir)
    , sA(m.sA)
    , sB(m.sB)
    , eA(m.sA + m.len - 1)
    , eB(m.sB + m.len - 1)
    , deltaApos(0)
    , Errors(0)
    , SimErrors(0)
    , NonAlphas(0)
  { }
  Alignment() = default;
  Alignment(const Alignment& rhs) = default;
  Alignment(Alignment&& rhs) = default;
  Alignment& operator=(const Alignment& rhs) = default;

  // Number of bases in alignment (in reference), counting deletions.
  long total() const {
    return std::abs(eA - sA) + 1 + std::count_if(delta.cbegin(), delta.cend(), [](long x) { return x < 0; });
  }

  inline double identity(const long t) const { return (double)(t - Errors) / t; }
  inline double identity() const { return identity(total()); }
  inline double similarity(const long t) const { return (double)(t - SimErrors) / t; }
  inline double similarity() const { return similarity(total()); }
  inline double stopity(const long t) const { return (double)NonAlphas / (2 * t); }
  inline double stopity() const { return stopity(total()); }

  struct stats_type {
    double identity;
    double similarity;
    double stopity;
  };
  stats_type stats() const {
    const long t = total();
    return { identity(t), similarity(t), stopity(t) };
  }
};

// Iterator over the error

// The type of error. INSERTION and DELETION apply to the
// reference. INSERTION means an there is an extra base in the
// reference., DELETION means there is one less base in the reference.
enum error_type { NONE, INSERTION, DELETION, MISMATCH };
struct error_description_type {
  error_type  type;             // The type of the error
  const char* ref;              // Pointer in the reference of the error
  const char* qry;              // Idem qry
  long        dst;              // Distance to previous indel or start of sequences

  // For a mismatch, *ref != *qry. For an INSERTION, qry is the
  // position in the qry sequence where to insert the extra base,
  // which is *ref. Idem for DELETION.
  //
  // dst is the distance to the previous indel. So dst >= 1 and dst-1
  // is equal to the number of matching and mismatching bases since the start or the
  // previous indel.
};
class error_iterator_type : public std::iterator<std::input_iterator_tag, error_description_type> {
  const Alignment&       m_al;
  error_description_type m_error;
  const char*            m_ref_end;
  size_t                 m_k;   // index in delta
public:
  // Create an iterator at beginning of error. ref and qry pointer
  // points char before the start of sequence (i.e. as in 1-base
  // indexing).
  error_iterator_type(const Alignment& al, const char* ref, const char* qry, size_t qry_len)
    : m_al(al)
    , m_error{ NONE, ref + al.sA - 1, qry + (al.dirB == 1 ? al.sB - 1 : (long)qry_len - al.sB + 2), 0 }
    , m_ref_end(ref + al.eA + 1)
    , m_k(0)
  { ++*this; }
  // Create an iterator at end
  error_iterator_type(const Alignment& al, const char* ref)
    : m_al(al)
    , m_error{NONE, ref + al.eA + 1, nullptr, 0}
    , m_ref_end(nullptr)
    , m_k(0)
  { }
  static char comp(char b) {
    switch(b) {
    case 'a': return 't'; case 'A': return 'T';
    case 'c': return 'g'; case 'C': return 'G';
    case 'g': return 'c'; case 'G': return 'C';
    case 't': return 'a'; case 'T': return 'A';
    default: return 'n';
    }
  }

  bool operator==(const error_iterator_type& rhs) const { return m_error.ref == rhs.m_error.ref; }
  bool operator!=(const error_iterator_type& rhs) const { return m_error.ref != rhs.m_error.ref; }
  const error_description_type& operator*() const { return m_error; }
  const error_description_type* operator->() const { return &m_error; }
  error_iterator_type& operator++();
  error_iterator_type operator++(int) {
    error_iterator_type res(*this);
    ++*this;
    return res;
  }
};

std::ostream& operator<<(std::ostream& os, const Alignment& al);

struct AscendingClusterSort
//-- For sorting clusters in ascending order of their sA coordinate
{
  bool operator() (const Cluster & pA, const Cluster & pB)
  {
    return ( pA.matches.front().sA < pB.matches.front().sA );
  }
};

struct merge_syntenys {
  const bool                     DO_DELTA;
  const bool                     DO_EXTEND;
  const bool                     TO_SEQEND;
  const bool                     DO_SHADOWS;
  const sw_align::aligner_buffer aligner;

  merge_syntenys(bool dd, bool de, bool ts, bool ds)
    : DO_DELTA(dd)
    , DO_EXTEND(de)
    , TO_SEQEND(ts)
    , DO_SHADOWS(ds)
    , aligner()
  { }

  merge_syntenys(bool dd, bool de, bool ts, bool ds, int break_len, int banding, int matrix_type)
    : DO_DELTA(dd)
    , DO_EXTEND(de)
    , TO_SEQEND(ts)
    , DO_SHADOWS(ds)
    , aligner(break_len, banding, matrix_type)
  { }

  // Process all syntenys in a container
  template<typename Container, typename FR2, typename ClustersOut, typename MatchesOut>
  void processSyntenys_each(Container& Syntenys, const FR2& Bf,
                            ClustersOut clusters, MatchesOut matches) const;
  template<typename Container, typename FR2, typename MatchesOut>
  void processSyntenys_each(Container& Syntenys, const FR2& Bf,
                            MatchesOut matches) const {
    processSyntenys_each(Syntenys, Bf, [](const Container& s, const FR2& Bf) { },
                         matches);
  }

  // Process all syntenys in a container in parallel
  template<typename Container, typename FR2, typename ClustersOut, typename MatchesOut>
  void processSyntenys_long_each(Container& Syntenys, const FR2& Bf,
                                 ClustersOut clusters, MatchesOut matches) const;
  // Process all syntenys in a container in parallel
  template<typename Container, typename FR2, typename MatchesOut>
  void processSyntenys_long_each(Container& Syntenys, const FR2& Bf,
                                 MatchesOut matches) const {
    processSyntenys_long_each(Syntenys, Bf, [](const Container& s, const FR2& Bf) { },
                              matches);
  }

  bool extendBackward(std::vector<Alignment> & Alignments, std::vector<Alignment>::iterator CurrAp,
                      std::vector<Alignment>::iterator TargetAp, const char * A, const char * B) const;

  void extendClusters(std::vector<Cluster> & Clusters,
                      const char* Aseq, const long Alen, const char* Bseq, const long Blen,
                      std::vector<Alignment>& Alignments) const;

  std::vector<Alignment> extendClusters(std::vector<Cluster> & Clusters,
                                        const char* Aseq, const long Alen, const char* Bseq, const long Blen) const {
    std::vector<Alignment> res;
    extendClusters(Clusters, Aseq, Alen, Bseq, Blen, res);
    return res;
  }

protected:
  bool extendForward(std::vector<Alignment>::iterator Ap, const char * A, long int targetA,
                     const char * B, long int targetB, unsigned int m_o) const;

  std::vector<Cluster>::iterator getForwardTargetCluster(std::vector<Cluster> & Clusters, std::vector<Cluster>::iterator CurrCp,
                                                         long int & targetA, long int & targetB) const;
  std::vector<Alignment>::iterator getReverseTargetAlignment(std::vector<Alignment> & Alignments,
                                                             std::vector<Alignment>::iterator CurrAp) const;
  void parseDelta(std::vector<Alignment> & Alignments,
                  const char* Aseq, const char* Bseq, const long Blen) const;
};

//-- Helper functions
inline void ignore_line(std::istream& is) {
  is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

//------------------------------------------------- Function Declarations ----//
bool Read_Sequence(std::istream& is, std::string& T, std::string& name);

void printDeltaAlignments(const std::vector<Alignment>& Alignments,
                          const std::string& AId, const long Alen,
                          const std::string& BId, const long Blen,
                          std::ostream& DeltaFile, const long minLen = 0);

template<typename FastaRecord>
inline void printDeltaAlignments(const std::vector<Alignment>& Alignments,
                          const FastaRecord& Af, const FastaRecord& Bf,
                          std::ostream& DeltaFile) {
  printDeltaAlignments(Alignments, Af.Id(), Af.len(), Bf.Id(), Bf.len(), DeltaFile);
}

// Print alignments in SAM format
template<typename FR1, typename FR2>
void printSAMAlignments(const std::vector<Alignment>& Alignments,
                        const FR1& A, const FR2& B,
                        std::ostream& SAMFile, bool long_format, const long minLen = 0);
std::string createCIGAR(const std::vector<long int>& ds, long int start, long int end, long int len, bool hard_clip = false);
std::string createMD(const Alignment& al, const char* ref,
                     const char* qry, size_t qry_len);

template<typename FastaRecord>
void printSyntenys(const std::vector<Synteny<FastaRecord> >& Syntenys, const FastaRecord& Bf, std::ostream& ClusterFile);

bool isShadowedCluster
(std::vector<Cluster>::const_iterator CurrCp,
 const std::vector<Alignment> & Alignments, std::vector<Alignment>::const_iterator Ap);

void __parseAbort
(const char *msg, const char* file, size_t line);
inline void __parseAbort(const std::string& s, const char* file, size_t line) {
  __parseAbort(s.c_str(), file, line);
}

#define parseAbort(msg) __parseAbort(msg, __FILE__, __LINE__);

inline long int revC
(long int Coord, long int Len)
  //  Reverse complement the given coordinate for the given length.

{
  assert (Len - Coord + 1 > 0);
  return (Len - Coord + 1);
}



//
// Implementation of templated methods
//
template<typename Container, typename FR2, typename ClustersOut, typename MatchesOut>
void merge_syntenys::processSyntenys_each(Container& Syntenys, const FR2& Bf,
                                          ClustersOut clusters, MatchesOut matches) const

//  For each syntenic region with clusters, extend the clusters to
//  expand total alignment coverage. Only should be called once all
//  the clusters for the contained syntenic regions have been stored
//  in the data structure. Frees the memory used by the the syntenic
//  regions once the output of extendClusters and flushSyntenys has
//  been produced.

{
  //-- For all the contained syntenys
  for(auto& CurrSp : Syntenys) {
    //-- If no clusters, ignore
    if(CurrSp.clusters.empty()) continue;
    //-- Extend clusters and create the alignment information
    //      alignments.clear();
    std::vector<Alignment> alignments;
    extendClusters (CurrSp.clusters, CurrSp.AfP->seq(), CurrSp.AfP->len(), Bf.seq(), Bf.len(), alignments);
    //-- Output the alignment data to the delta file
    matches(std::move(alignments), *CurrSp.AfP, Bf);
  }

  //-- Create the cluster information
  clusters(Syntenys, Bf);
  Syntenys.clear();
}

template<typename Container, typename FR2, typename ClustersOut, typename MatchesOut>
void merge_syntenys::processSyntenys_long_each(Container& Syntenys, const FR2& Bf,
                                               ClustersOut clusters, MatchesOut matches) const

//  For each syntenic region with clusters, extend the clusters to
//  expand total alignment coverage. Only should be called once all
//  the clusters for the contained syntenic regions have been stored
//  in the data structure. Frees the memory used by the the syntenic
//  regions once the output of extendClusters and flushSyntenys has
//  been produced.

{
  std::atomic<size_t> pos(0);
  const auto end = Syntenys.end();

//#pragma omp parallel
  {
    auto CurrSp = Syntenys.begin();
    size_t prev = 0;

    while(true) {
      const size_t npos = pos++;
      for( ; prev < npos && CurrSp != end; ++prev, ++CurrSp) ; // Advance
      if(CurrSp == end) break;
      if(CurrSp->clusters.empty()) continue;
      //-- Extend clusters and create the alignment information
      //      alignments.clear();
      std::vector<Alignment> alignments;
      extendClusters (CurrSp->clusters, CurrSp->AfP->seq(), CurrSp->AfP->len(), Bf.seq(), Bf.len(), alignments);
      //-- Output the alignment data to the delta file
      matches(std::move(alignments), *CurrSp->AfP, Bf);
    }
  }

  //-- Create the cluster information
  clusters(Syntenys, Bf);
  Syntenys.clear();
}

template<typename FastaRecord>
void printSyntenys(const std::vector<Synteny<FastaRecord> > & Syntenys, const FastaRecord& Bf, std::ostream& ClusterFile)

//  Simply output the synteny/cluster information generated by the mgaps
//  program. However, now the coordinates reference their appropriate
//  reference sequence, and the reference sequecne header is added to
//  the appropriate lines. Free the memory used by Syntenys once the
//  data is successfully output to the file.

{
  if ( ClusterFile ) {
    for(const auto& Sp : Syntenys) { // each syntenys
      ClusterFile << '>' << Sp.AfP->Id() << ' ' << Bf.Id() << ' '
                  << Sp.AfP->len() << ' ' << Bf.len() << '\n';

      for (const auto& Cp : Sp.clusters) { // each clusters
        ClusterFile << std::setw(2) << FORWARD_CHAR << ' ' << std::setw(2) << Cp.dirB << '\n';

        for (auto Mp = Cp.matches.cbegin( ); Mp != Cp.matches.cend( ); ++Mp ) { // each match
            ClusterFile << std::setw(8) << Mp->sA << ' '
                        << std::setw(8) << (Cp.dirB == FORWARD_CHAR ? Mp->sB : revC(Mp->sB, Bf.len())) << ' '
                        << std::setw(6) << Mp->len;
          if ( Mp != Cp.matches.cbegin( ) )
            ClusterFile << std::setw(6) << (Mp->sA - (Mp - 1)->sA - (Mp - 1)->len) << ' '
                        << std::setw(6) << (Mp->sB - (Mp - 1)->sB - (Mp - 1)->len) << '\n';
          else
            ClusterFile << "     -      -\n";
        }
      }
    }
  }
}

template<typename FR1, typename FR2>
void printSAMAlignments(const std::vector<Alignment>& Alignments,
                        const FR1& A, const FR2& B,
                        std::ostream& SAMFile, bool long_format,
                        const long minLen) {
  const char* mapq = Alignments.size() > 1 ? "\t10\t" : "\t30\t";
  bool hard_clip = false;
  for(const auto& Al : Alignments) {
    if(std::abs(Al.eA - Al.sA) < minLen && std::abs(Al.eB - Al.sB) < minLen)
      continue;
    const bool fwd = Al.dirB == FORWARD_CHAR;
    SAMFile << B.Id() << '\t'
            << ((hard_clip ? 0x800 : 0) | (fwd ? 0 : 0x10)) << '\t'
            << A.Id() << '\t' << Al.sA
            << mapq
            << createCIGAR(Al.delta, Al.sB, Al.eB, B.len(), hard_clip)
            << "\t*\t0\t0\t";
    if(long_format) {
      const auto start = hard_clip ? Al.sB : 1;
      const auto len   = hard_clip ? Al.eB - start + 1 : B.len();
      SAMFile.write(B.seq() + start, len);
    } else {
      SAMFile << '*';
    }
    SAMFile << "\t*\tNM:i:" << Al.Errors;
    if(long_format)
      SAMFile << "\tMD:Z:" << createMD(Al, A.seq(), B.seq(), B.len());
    SAMFile << '\n';
    hard_clip = true;
  }
}


} // namespace postnuc
} // namespace mummer

#endif /* __MUMMER_POSTNUC_H__ */
