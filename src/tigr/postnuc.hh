#ifndef __MUMMER_POSTNUC_H__
#define __MUMMER_POSTNUC_H__

#include <string>
#include <vector>
#include <iostream>
#include <limits>


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
};


struct Synteny
//-- An ordered list of clusters between two sequences A and B
{
  FastaRecord *   AfP;          // a pointer to the reference sequence record
  FastaRecord     Bf;           // the query sequence record (w/o the sequence)
  std::vector<Cluster> clusters;     // the ordered set of clusters between A and B
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

  void processSyntenys(std::vector<Synteny> & Syntenys,
                       FastaRecord * Af,
                       std::istream& QryFile, std::ostream& ClusterFile, std::ostream& DeltaFile);
  bool extendBackward(std::vector<Alignment> & Alignments, std::vector<Alignment>::iterator CurrAp,
                      std::vector<Alignment>::iterator TargetAp, const char * A, const char * B);
  void extendClusters(std::vector<Cluster> & Clusters,
                      const FastaRecord * Af, const FastaRecord * Bf, std::ostream& DeltaFile);
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

void addNewAlignment
(std::vector<Alignment> & Alignments, std::vector<Cluster>::iterator Cp,
 std::vector<Match>::iterator Mp);

bool extendForward
(std::vector<Alignment>::iterator Ap, const char * A, long int targetA,
 const char * B, long int targetB, unsigned int m_o);

void flushAlignments
(std::vector<Alignment> & Alignments,
 const FastaRecord * Af, const FastaRecord * Bf,
 std::ostream& DeltaFile);

void flushSyntenys
(std::vector<Synteny> & Syntenys, std::ostream& ClusterFile);

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


} // namespace postnuc
} // namespace mummer

#endif /* __MUMMER_POSTNUC_H__ */
