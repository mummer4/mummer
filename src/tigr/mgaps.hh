#ifndef __MUMMER_MGAPS_H__
#define __MUMMER_MGAPS_H__

namespace mummer {
namespace mgaps {
// Represent a match
struct  Match_t {
  long int     Start1, Start2, Len;
  long int     Simple_Score;
  long int     Simple_From;
  long int     Simple_Adj;
  unsigned int cluster_id:30;
  unsigned int Good:1;
  unsigned int Tentative:1;
  Match_t() = default;
  Match_t(long int S1, long int S2, long int L) : Start1(S1), Start2(S2), Len(L), Good(true), Tentative(false) { }
};


// Union find data structures. Indices MUST be in the range [1, s],
// where s is the size given to reset (no checks).
struct UnionFind {
  std::vector<int> m_UF;

  void reset(size_t s);         // Reset to given size
  int  find(int a);             //  Return the id of the set containing  a  in  UF .
  void union_sets(int a, int b); //  Union the sets whose id's are  a  and  b  in  UF .
};

struct ClusterMatches {
  const int      Fixed_Separation;
  const long int Max_Separation;
  const long int Min_Output_Score;
  const double   Separation_Factor;
  const bool     Use_Extents;

  ClusterMatches(int fs, long int ms, long int mo, double sf, bool ue)
    : Fixed_Separation(fs)
    , Max_Separation(ms)
    , Min_Output_Score(mo)
    , Separation_Factor(sf)
    , Use_Extents(ue)
  { }

  typedef std::vector<std::vector<Match_t>> clusters_type;

  int  Process_Matches(Match_t * A, UnionFind& UF, int N, clusters_type& clusters);
  int  Process_Cluster(Match_t * A, int N, clusters_type& cluster);
  void Print_Cluster(const clusters_type& clusters, const char* label);
};


} // namespace mgaps
} // namespace mummer

#endif /* __MUMMER_MGAPS_H__ */
