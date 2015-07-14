#ifndef __MUMMER_POSTNUC_H__
#define __MUMMER_POSTNUC_H__

#include <string>
#include <vector>
#include <iostream>
#include <limits>
#include <cstring>
#include <memory>
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


class FastaRecord
//-- The essential data of a sequence
{
  std::string m_Id;               // the fasta ID header tag
  long int    m_len;              // the length of the sequence
  std::string m_seq;              // the sequence data

public:
  FastaRecord() = default;
  const std::string& Id() const { return m_Id; }
  long int len() const { return m_len; }
  const char* seq() const { return m_seq.c_str(); }

  long int& len_w() { return m_len; }
  std::string& Id_w() { return m_Id; }

  friend bool Read_Sequence(std::istream& is, FastaRecord& record);
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

  template<typename ClustersOut, typename MatchesOut>
  void processSyntenys_each(std::vector<Synteny>& Syntenys, const FastaRecord& Bf,
                            ClustersOut clusters, MatchesOut matches);
  bool extendBackward(std::vector<Alignment> & Alignments, std::vector<Alignment>::iterator CurrAp,
                      std::vector<Alignment>::iterator TargetAp, const char * A, const char * B);

  template<typename MatchesOut>
  void extendClusters(std::vector<Cluster> & Clusters,
                      const FastaRecord * Af, const FastaRecord * Bf, MatchesOut matches);
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
bool Read_Sequence(std::istream& is, FastaRecord& record);

void printDeltaAlignments(const std::vector<Alignment> & Alignments,
                          const FastaRecord * Af, const FastaRecord * Bf,
                          std::ostream& DeltaFile);

void printSyntenys(const std::vector<Synteny>& Syntenys, const FastaRecord& Bf, std::ostream& ClusterFile);

bool extendForward
(std::vector<Alignment>::iterator Ap, const char * A, long int targetA,
 const char * B, long int targetB, unsigned int m_o);

void flushAlignments
(std::vector<Alignment> & Alignments,
 const FastaRecord * Af, const FastaRecord * Bf,
 std::ostream& DeltaFile);

void flushSyntenys
(std::vector<Synteny> & Syntenys, const FastaRecord& Bf, std::ostream& ClusterFile);

std::vector<Cluster>::iterator getForwardTargetCluster
(std::vector<Cluster> & Clusters, std::vector<Cluster>::iterator CurrCp,
 long int & targetA, long int & targetB);

std::vector<Alignment>::iterator getReverseTargetAlignment
(std::vector<Alignment> & Alignments, std::vector<Alignment>::iterator CurrAp);

bool isShadowedCluster
(std::vector<Cluster>::iterator CurrCp,
 std::vector<Alignment> & Alignments, std::vector<Alignment>::iterator Ap);

void __parseAbort
(const char *msg, const char* file, size_t line);
inline void __parseAbort(const std::string& s, const char* file, size_t line) {
  __parseAbort(s.c_str(), file, line);
}

#define parseAbort(msg) __parseAbort(msg, __FILE__, __LINE__);

void parseDelta
(std::vector<Alignment> & Alignments,
 const FastaRecord * Af, const FastaRecord *Bf);

inline long int revC
(long int Coord, long int Len);


//
// Implementation of templated methods
//
template<typename ClustersOut, typename MatchesOut>
void merge_syntenys::processSyntenys_each(std::vector<Synteny>& Syntenys, const FastaRecord& Bf,
                                          ClustersOut clusters, MatchesOut matches)

//  For each syntenic region with clusters, extend the clusters to
//  expand total alignment coverage. Only should be called once all
//  the clusters for the contained syntenic regions have been stored
//  in the data structure. Frees the memory used by the the syntenic
//  regions once the output of extendClusters and flushSyntenys has
//  been produced.

{
  //-- For all the contained syntenys
  for(auto CurrSp : Syntenys) {
      //-- If no clusters, ignore
      if(CurrSp.clusters.empty()) continue;
      //-- Extend clusters and create the alignment information
      extendClusters (CurrSp.clusters, CurrSp.AfP, &Bf, matches);
  }

  //-- Create the cluster information
  clusters(Syntenys, Bf);
  Syntenys.clear();
}

template<typename MatchesOut>
void merge_syntenys::extendClusters(std::vector<Cluster> & Clusters,
                                    const FastaRecord * Af, const FastaRecord * Bf, MatchesOut matches)

//  Connect all the matches in every cluster between sequences A and B.
//  Also, extend alignments off of the front and back of each cluster to
//  expand total alignment coverage. When these extensions encounter an
//  adjacent cluster, fuse the two regions to create one single
//  encompassing region. This routine will create alignment objects from
//  these extensions and output the resulting delta information to the
//  delta output file.

{
  //-- Sort the clusters (ascending) by their start coordinate in sequence A
  sort (Clusters.begin( ), Clusters.end( ), AscendingClusterSort( ));


  //-- If no delta file is requested
  if ( ! DO_DELTA )
    return;


  bool target_reached = false;         // reached the adjacent match or cluster

  const char * A, * B;                 // the sequences A and B
  std::unique_ptr<char[]> Brev;          // the reverse complement of B

  unsigned int m_o;
  long int targetA, targetB;           // alignment extension targets in A and B

  std::vector<Match>::iterator Mp;          // match pointer

  std::vector<Cluster>::iterator PrevCp;    // where the extensions last left off
  std::vector<Cluster>::iterator CurrCp;    // the current cluster being extended
  std::vector<Cluster>::iterator TargetCp = Clusters.end( );  // the target cluster

  std::vector<Alignment> Alignments;        // the vector of alignment objects
  std::vector<Alignment>::iterator CurrAp = Alignments.begin( );   // current align
  std::vector<Alignment>::iterator TargetAp;                // target align


  //-- Extend each cluster
  A = Af->seq();
  PrevCp = Clusters.begin( );
  CurrCp = Clusters.begin( );
  while ( CurrCp < Clusters.end( ) ) {
      if ( DO_EXTEND ) {
        if ( ! target_reached ) //-- Ignore if shadowed or already extended
          if ( CurrCp->wasFused ||
               (!DO_SHADOWS && isShadowedCluster (CurrCp, Alignments, CurrAp)) ) {
            CurrCp->wasFused = true;
            CurrCp = ++ PrevCp;
            continue;
          }
      }

      //-- Pick the right directional sequence for B
      if ( CurrCp->dirB == FORWARD_CHAR )
        B = Bf->seq();
      else if ( Brev )
        B = Brev.get();
      else {
        Brev.reset(new char[Bf->len() + 2]);
        strcpy ( Brev.get() + 1, Bf->seq() + 1 );
        Brev[0] = '\0';
        Reverse_Complement (Brev.get(), 1, Bf->len());
        B = Brev.get();
      }

      //-- Extend each match in the cluster
      for ( Mp = CurrCp->matches.begin( ); Mp < CurrCp->matches.end( ); Mp ++ ) {
        //-- Try to extend the current match backwards
        if ( target_reached ) {
          //-- Merge with the previous match
          if ( CurrAp->eA != Mp->sA  ||  CurrAp->eB != Mp->sB ) {
            if ( Mp >= CurrCp->matches.end( ) - 1 ) {
              std::cerr << "ERROR: Target match does not exist, please\n"
                   << "       file a bug report\n";
              exit (EXIT_FAILURE);
            }
            continue;
          }
          CurrAp->eA += Mp->len - 1;
          CurrAp->eB += Mp->len - 1;
        } else { //-- Create a new alignment object
          Alignments.push_back({ *Mp, CurrCp->dirB } );
          CurrAp = Alignments.end( ) - 1;

          if ( DO_EXTEND  ||  Mp != CurrCp->matches.begin ( ) ) {
            //-- Target the closest/best alignment object
            TargetAp = getReverseTargetAlignment (Alignments, CurrAp);

            //-- Extend the new alignment object backwards
            if ( extendBackward (Alignments, CurrAp, TargetAp, A, B) )
              CurrAp = TargetAp;
          }
        }

          m_o = FORWARD_ALIGN;

          //-- Try to extend the current match forwards
          if ( Mp < CurrCp->matches.end( ) - 1 ) {
            //-- Target the next match in the cluster
            targetA = (Mp + 1)->sA;
            targetB = (Mp + 1)->sB;

            //-- Extend the current alignment object forward
            target_reached = extendForward (CurrAp, A, targetA, B, targetB, m_o);
          } else if ( DO_EXTEND ) {
            targetA = Af->len();
            targetB = Bf->len();

            //-- Target the closest/best match in a future cluster
            TargetCp = getForwardTargetCluster (Clusters, CurrCp, targetA, targetB);
            if ( TargetCp == Clusters.end( ) ) {
              m_o |= OPTIMAL_BIT;
              if ( TO_SEQEND )
                m_o |= SEQEND_BIT;
            }

            //-- Extend the current alignment object forward
            target_reached = extendForward (CurrAp, A, targetA, B, targetB, m_o);
          }
      }
      if ( TargetCp == Clusters.end( ) )
        target_reached = false;

      CurrCp->wasFused = true;

      if ( target_reached == false )
        CurrCp = ++ PrevCp;
      else
        CurrCp = TargetCp;
    }

#ifdef _DEBUG_ASSERT
  validateData (Alignments, Clusters, Af, Bf);
#endif

  //-- Generate the error counts
  parseDelta(Alignments, Af, Bf);
  //-- Output the alignment data to the delta file
  matches(Alignments, Af, Bf);
  Alignments.clear( );
}


} // namespace postnuc
} // namespace mummer

#endif /* __MUMMER_POSTNUC_H__ */
