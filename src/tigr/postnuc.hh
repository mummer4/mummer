#ifndef __MUMMER_POSTNUC_H__
#define __MUMMER_POSTNUC_H__

#include <string>
#include <vector>
#include <iostream>
#include <limits>
#include <cstring>
#include <memory>
#include <iomanip>
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
  { }
  Alignment(Alignment&& rhs) = default;
};


struct AscendingClusterSort
//-- For sorting clusters in ascending order of their sA coordinate
{
  bool operator() (const Cluster & pA, const Cluster & pB)
  {
    return ( pA.matches.front().sA < pB.matches.front().sA );
  }
};

struct merge_syntenys {
  const bool DO_DELTA;
  const bool DO_EXTEND;
  const bool TO_SEQEND;
  const bool DO_SHADOWS;

  merge_syntenys(bool dd, bool de, bool ts, bool ds)
    : DO_DELTA(dd)
    , DO_EXTEND(de)
    , TO_SEQEND(ts)
    , DO_SHADOWS(ds)
  { }

  template<typename FastaRecord, typename ClustersOut, typename MatchesOut>
  void processSyntenys_each(std::vector<Synteny<FastaRecord> >& Syntenys, const FastaRecord& Bf,
                            ClustersOut clusters, MatchesOut matches);
  bool extendBackward(std::vector<Alignment> & Alignments, std::vector<Alignment>::iterator CurrAp,
                      std::vector<Alignment>::iterator TargetAp, const char * A, const char * B);

  void extendClusters(std::vector<Cluster> & Clusters,
                      const char* Aseq, const long Alen, const char* Bseq, const long Blen,
                      std::vector<Alignment>& Alignments);

  std::vector<Alignment> extendClusters(std::vector<Cluster> & Clusters,
                                        const char* Aseq, const long Alen, const char* Bseq, const long Blen) {
    std::vector<Alignment> res;
    extendClusters(Clusters, Aseq, Alen, Bseq, Blen, res);
    return res;
  }
};

//-- Helper functions
template<typename stream_type>
void open_stream(stream_type& st, const std::string& name) {
  st.open(name);
  if(!st.good()) {
    std::cerr << "ERROR: Could not open file " << name << std::endl;
    exit(EXIT_FAILURE);
  }
}

inline void open_ofstream(std::ofstream& os, const std::string& name) {
  open_stream<std::ofstream>(os, name);
}

inline void open_ifstream(std::ifstream& os, const std::string& name) {
  open_stream<std::ifstream>(os, name);
}

inline void ignore_line(std::istream& is) {
  is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

//------------------------------------------------- Function Declarations ----//
bool Read_Sequence(std::istream& is, std::string& T, std::string& name);

template<typename FastaRecord>
void printDeltaAlignments(const std::vector<Alignment>& Alignments,
                          const FastaRecord& Af, const FastaRecord& Bf,
                          std::ostream& DeltaFile);

template<typename FastaRecord>
void printSyntenys(const std::vector<Synteny<FastaRecord> >& Syntenys, const FastaRecord& Bf, std::ostream& ClusterFile);

bool extendForward
(std::vector<Alignment>::iterator Ap, const char * A, long int targetA,
 const char * B, long int targetB, unsigned int m_o);

std::vector<Cluster>::iterator getForwardTargetCluster
(std::vector<Cluster> & Clusters, std::vector<Cluster>::iterator CurrCp,
 long int & targetA, long int & targetB);

std::vector<Alignment>::iterator getReverseTargetAlignment
(std::vector<Alignment> & Alignments, std::vector<Alignment>::iterator CurrAp);

bool isShadowedCluster
(std::vector<Cluster>::const_iterator CurrCp,
 const std::vector<Alignment> & Alignments, std::vector<Alignment>::const_iterator Ap);

void __parseAbort
(const char *msg, const char* file, size_t line);
inline void __parseAbort(const std::string& s, const char* file, size_t line) {
  __parseAbort(s.c_str(), file, line);
}

#define parseAbort(msg) __parseAbort(msg, __FILE__, __LINE__);

void parseDelta(std::vector<Alignment> & Alignments,
                const char* Aseq, const char* Bseq, const long Blen);

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
template<typename FastaRecord, typename ClustersOut, typename MatchesOut>
void merge_syntenys::processSyntenys_each(std::vector<Synteny<FastaRecord> >& Syntenys, const FastaRecord& Bf,
                                          ClustersOut clusters, MatchesOut matches)

//  For each syntenic region with clusters, extend the clusters to
//  expand total alignment coverage. Only should be called once all
//  the clusters for the contained syntenic regions have been stored
//  in the data structure. Frees the memory used by the the syntenic
//  regions once the output of extendClusters and flushSyntenys has
//  been produced.

{
  //  std::vector<Alignment> alignments;

  //-- For all the contained syntenys
  for(auto CurrSp : Syntenys) {
      //-- If no clusters, ignore
      if(CurrSp.clusters.empty()) continue;
      //-- Extend clusters and create the alignment information
      //      alignments.clear();
      std::vector<Alignment> alignments;
      extendClusters (CurrSp.clusters, CurrSp.AfP->seq(), CurrSp.AfP->len(), Bf.seq(), Bf.len(), alignments);
      //-- Output the alignment data to the delta file
      matches(alignments, *CurrSp.AfP, Bf);
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

template<typename FastaRecord>
void printDeltaAlignments(const std::vector<Alignment> & Alignments,
                          const FastaRecord& Af, const FastaRecord& Bf,
                          std::ostream& DeltaFile)

//  Simply output the delta information stored in Alignments to the
//  given delta file. Free the memory used by Alignments once the
//  data is successfully output to the file.

{
  DeltaFile << '>' << Af.Id() << ' ' << Bf.Id() << ' ' << Af.len() << ' ' << Bf.len() << '\n';

  for(const auto& A : Alignments) {
    const bool fwd = A.dirB == FORWARD_CHAR;
    DeltaFile << A.sA << ' ' << A.eA << ' '
              << (fwd ? A.sB : revC(A.sB, Bf.len())) << ' '
              << (fwd ? A.eB : revC(A.eB, Bf.len())) << ' '
              << A.Errors << ' ' << A.SimErrors << ' ' << A.NonAlphas
              << '\n';

    for(const auto& D : A.delta)
      DeltaFile << D << '\n';
    DeltaFile << "0\n";
  }
}

} // namespace postnuc
} // namespace mummer

#endif /* __MUMMER_POSTNUC_H__ */
