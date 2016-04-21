#ifndef __MUMMER_MGAPS_H__
#define __MUMMER_MGAPS_H__

#include <iostream>
#include <cassert>
#include <mummer/dset.hpp>
#include <mummer/openmp_qsort.hpp>

namespace mummer {
namespace mgaps {
// Represent a match
struct  Match_t {
  long int     Start1, Start2, Len; // Start1 and Start2 are 1-based.
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
typedef std::vector<Match_t>      cluster_type;
typedef std::vector<cluster_type> clusters_type;

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

  template<typename Output>
  int Cluster_each(Match_t * A, UnionFind& UF, int N, Output out) const;

  // Like Cluster_each, but adapted for long query with many matches.
  template<typename Output>
  int Cluster_each_long(Match_t * A, int N, Output out) const;

  //  Process matches  A [1 .. N]  and append them to clusters
  int  Process_Matches(Match_t * A, UnionFind& UF, int N, clusters_type& clusters) const {
    return Cluster_each(A, UF, N, [&](cluster_type&& cl) { clusters.push_back(std::move(cl)); });
  }

  static void Print_Cluster(const cluster_type& cluster, const char* label, std::ostream& os = std::cout);
  static void Print_Clusters(const clusters_type& clusters, const char* label, std::ostream& os = std::cout);


protected:
  template<typename Output>
  int  Process_Cluster(Match_t * A, int N, Output out) const;

  //  Remove from  A [0 .. (N - 1)]  any matches that are internal to a repeat,
  static int Filter_Matches(Match_t* A, const int N);

  // Matches ordering
  static inline bool By_Start2(const Match_t& A, const Match_t& B) {
    return (A.Start2 < B.Start2) || (A.Start2 == B.Start2 && A.Start1 < B.Start1);
  }

  static bool By_Cluster(const Match_t& A, const Match_t& B) {
    return (A.cluster_id < B.cluster_id) ||
      (A.cluster_id == B.cluster_id && By_Start2(A, B));
  }
};

//
// Implementation of templated methods
//

template<typename Output>
int ClusterMatches::Cluster_each(Match_t * A, UnionFind& UF, int N, Output out) const {
  //  Process matches  A [1 .. N]  and output them after
  //  a line containing  label .

  //  Use Union-Find to create connected-components based on
  //  separation and similar diagonals between matches
  UF.reset(N);

  std::sort(A + 1, A + N + 1, By_Start2);
  N = Filter_Matches (A + 1, N);

  for  (int i = 1;  i < N;  i ++) {
    long int i_end  = A [i] . Start2 + A [i] . Len;
    long int i_diag = A [i] . Start2 - A [i] . Start1;

    for  (int j = i + 1;  j <= N;  j ++) {
      long int sep = A [j] . Start2 - i_end;
      if  (sep > Max_Separation)
        break;

      long int diag_diff = std::abs ((A [j] . Start2 - A [j] . Start1) - i_diag);
      if  (diag_diff <= std::max(Fixed_Separation, (int)(Separation_Factor * sep)))
        UF.union_sets(UF.find(i), UF.find(j));
    }
  }

  //  Set the cluster id of each match and reset Good flag
  for  (int i = 1;  i <= N;  i ++) {
    A [i] . cluster_id = UF.find (i);
    assert(A[i].cluster_id > 0);
    A[i].Good = false;
  }
  std::sort(A + 1, A + N + 1, By_Cluster);

  // Determine and process clusters
  int cluster_size, print_ct = 0;
  for (int i = 1;  i <= N;  i += cluster_size) {
    int j;
    for  (j = i + 1;  j <= N && A [i] . cluster_id == A [j] . cluster_id;  j ++)
      ;
    cluster_size = j - i;
    print_ct += Process_Cluster (A + i, cluster_size, out);
  }
  return print_ct;
}

template<typename Output>
int ClusterMatches::Cluster_each_long(Match_t * A, int N, Output out) const {
  //  Process matches  A [1 .. N]  and output them after
  //  a line containing  label .

  //  Use Union-Find to create connected-components based on
  //  separation and similar diagonals between matches
  DisjointSets UF(N + 1);

  openmp_qsort(A + 1, A + N + 1, By_Start2);
  N = Filter_Matches (A + 1, N);

#pragma omp parallel for schedule(dynamic)
  for  (int i = 1;  i < N;  i ++) {
    long int i_end  = A [i] . Start2 + A [i] . Len;
    long int i_diag = A [i] . Start2 - A [i] . Start1;

    for  (int j = i + 1;  j <= N;  j ++) {
      long int sep = A [j] . Start2 - i_end;
      if  (sep > Max_Separation)
        break;

      long int diag_diff = std::abs ((A [j] . Start2 - A [j] . Start1) - i_diag);
      if  (diag_diff <= std::max(Fixed_Separation, (int)(Separation_Factor * sep)))
        UF.union_sets(UF.find(i), UF.find(j));
    }
  }

  //  Set the cluster id of each match and reset Good flag
#pragma omp parallel for
  for  (int i = 1;  i <= N;  i ++) {
    A [i] . cluster_id = UF.find (i);
    assert(A[i].cluster_id > 0);
    A[i].Good = false;
  }
  openmp_qsort(A + 1, A + N + 1, By_Cluster);

  // Determine and process clusters
  int cluster_size, print_ct = 0;
#pragma omp parallel
  {
#pragma omp single
    for (int i = 1;  i <= N;  i += cluster_size) {
      int j;
      for  (j = i + 1;  j <= N && A [i] . cluster_id == A [j] . cluster_id;  j ++)
        ;
      cluster_size = j - i;
#pragma omp task firstprivate(i, cluster_size) shared(out)
      print_ct += Process_Cluster (A + i, cluster_size, out);
    }
  }

  return print_ct;
}

template<typename Output>
int ClusterMatches::Process_Cluster(Match_t * A, int N, Output out) const {
//  Process the cluster of matches in  A [0 .. (N - 1)]  and output them
//  after a line containing  label .  Return the number of clusters
//  printed.
  int  count = 0;

  while(N > 0) {
    std::vector<Match_t> cluster; // Potential cluster

    for  (int i = 0;  i < N;  i ++) {
      A [i] . Simple_Score = A [i] . Len;
      A [i] . Simple_Adj = 0;
      A [i] . Simple_From = -1;
      for  (int j = 0;  j < i;  j ++) {
        const long int Olap1 = A [j] . Start1 + A [j] . Len - A [i] . Start1;
        const long int Olap2 = A [j] . Start2 + A [j] . Len - A [i] . Start2;
        const long int Olap = std::max(std::max((long)0, Olap1), Olap2);

          // penalize off diagonal matches
        const long int Pen = Olap + std::abs ( (A [i] . Start2 - A [i] . Start1) -
                                               (A [j] . Start2 - A [j] . Start1) );

        if  (A [j] . Simple_Score + A [i] . Len - Pen > A [i] . Simple_Score) {
          A [i] . Simple_From = j;
          A [i] . Simple_Score = A [j] . Simple_Score + A [i] . Len - Pen;
          A [i] . Simple_Adj = Olap;
        }
      }
    }

    int best = 0;
    for  (int i = 1;  i < N;  i ++)
      if  (A [i] . Simple_Score > A [best] . Simple_Score)
        best = i;
    long int total = 0;
    long int hi    = LONG_MIN;
    long int lo    = LONG_MAX;
    for  (int i = best;  i >= 0;  i = A [i] . Simple_From) {
      A [i] . Good = true;
      total += A [i] . Len;
      hi = std::max(hi, A[i].Start1 + A[i].Len);
      lo = std::min(lo, A[i].Start1);
    }
    const long int score = Use_Extents ? hi - lo : total;

    if  (score >= Min_Output_Score) {
      count ++;
      for  (int i = 0;  i < N;  i ++)
        if  (A [i] . Good)
          cluster.push_back(A[i]);
      out(std::move(cluster));
    }

    // Compact match array
    const auto new_end = std::remove_if(A, A + N, [](const Match_t& m) { return m.Good; });
    N = new_end - A;
  }

  return count;
}


} // namespace mgaps
} // namespace mummer

#endif /* __MUMMER_MGAPS_H__ */
