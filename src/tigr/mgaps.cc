/* Programmer:  A. Delcher
*        File:  mgaps.cc
*
*  This program reads lists of unique matches between a sequence of strings
*  and a reference string.  For each string in the sequence, it clusters
*  the matches together into groups that may represent longer, inexact
*  matches.
*/

#include <cassert>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

#include  "tigrinc.hh"
#include "mgaps.hh"

namespace mummer {
namespace mgaps {


  // If true use end minus start as length of cluster instead of
  // sum of component lengths

// UnionFind
void UnionFind::reset(size_t s) {
  m_UF.resize(0);
  m_UF.resize(s + 1, -1);
  assert(m_UF[0] == -1 && m_UF[s] == -1);
}

int UnionFind::find(int a) { //  Return the id of the set containing  a  in  UF .
  if(m_UF [a] < 0) return  a;
  int i;
  for(i = a;  m_UF [i] > 0;  i = m_UF [i]) ;
  for(int k, j = a;  m_UF [j] != i;  j = k) {
    k = m_UF [j];
    m_UF [j] = i;
  }
  return  i;
}

void UnionFind::union_sets(int a, int b) { //  Union the sets whose id's are  a  and  b  in  UF .
  if(a == b) return;
  assert (m_UF [a] < 0 && m_UF [b] < 0);

  if(m_UF [a] < m_UF [b]) {
    m_UF [a] += m_UF [b];
    m_UF [b]  = a;
  } else {
    m_UF [b] += m_UF [a];
    m_UF [a]  = b;
  }
}

// Matches ordering
static bool By_Start2(const Match_t& A, const Match_t& B) {
  return (A.Start2 < B.Start2) || (A.Start2 == B.Start2 && A.Start1 < B.Start1);
}

static bool By_Cluster(const Match_t& A, const Match_t& B) {
  return (A.cluster_id < B.cluster_id) ||
    (A.cluster_id == B.cluster_id && By_Start2(A, B));
}


static void  Filter_Matches(Match_t * A, int & N) {
//  Remove from  A [0 .. (N - 1)]  any matches that are internal to a repeat,
//  e.g., if seq1 has 27 As and seq2 has 20 then the first and
//  last matches will be kept, but the 6 matches in the middle will
//  be eliminated.  Also combine overlapping matches on the same
//  diagonal.  Pack all remaining matches into the front of  A  and
//  reduce the value of  N  if any matches are removed.
//  Matches in  A  *MUST* be sorted by  Start2  value.
  for  (int i = 0;  i < N - 1;  i ++) {
    if  (! A[i].Good) continue;

    const int i_diag = A[i].Start2 - A[i].Start1;
    int       i_end  = A[i].Start2 + A[i].Len;

    for  (int j = i + 1;  j < N && A[j].Start2 <= i_end;  j ++) {
      assert (A[i].Start2 <= A[j].Start2);
      if  (! A[j].Good) continue;
      int j_diag = A[j].Start2 - A[j].Start1;
      if  (i_diag == j_diag) {
        int  j_extent = A[j].Len + A[j].Start2 - A[i].Start2;
        if  (j_extent > A[i].Len) {
          A[i].Len = j_extent;
          i_end = A[i].Start2 + j_extent;
        }
        A[j].Good = false;
      } else if  (A[i].Start1 == A[j].Start1) {
        int olap = A[i].Start2 + A[i].Len - A[j].Start2;
        if  (A[i].Len < A[j].Len) {
          if  (olap >=  A[i].Len / 2) {
            A[i].Good = false;
            break;
          }
        } else if  (A[j].Len < A[i].Len) {
          if  (olap >= A[j].Len / 2)
            A[j].Good = false;
        } else {
          if  (olap >= A[i].Len / 2) {
            A[j].Tentative = true;
            if  (A[i].Tentative) {
              A[i].Good = false;
              break;
            }
          }
        }
      } else if  (A[i].Start2 == A[j].Start2) {
        int olap = A[i].Start1 + A[i].Len - A[j].Start1;
        if  (A[i].Len < A[j].Len) {
          if  (olap >=  A[i].Len / 2) {
            A[i].Good = false;
            break;
          }
        } else if  (A[j].Len < A[i].Len) {
          if  (olap >= A[j].Len / 2)
            A[j].Good = false;

        } else {
          if  (olap >= A[i].Len / 2) {
            A[j].Tentative = true;
            if  (A[i].Tentative) {
              A[i].Good = false;
              break;
            }
          }
        }
      }
    }
  }

  int j = 0;
  for  (int i = 0;  i < N;  i ++) {
    if  (A[i].Good) {
      if  (i != j)
        A[j] = A[i];
      j ++;
    }
  }
  N = j;
}

int ClusterMatches::Process_Cluster(Match_t * A, int N, std::vector<std::vector<Match_t>>& out) {
//  Process the cluster of matches in  A [0 .. (N - 1)]  and output them
//  after a line containing  label .  Return the number of clusters
//  printed.
  int  count = 0;

  do {
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
      std::vector<Match_t>* cluster = nullptr;
      for  (int i = 0;  i < N;  i ++)
        if  (A [i] . Good) {
          if(!cluster) {
            out.push_back(std::vector<Match_t>());
            cluster = &out.back();
          }
          assert(cluster != nullptr);
          cluster->push_back(A[i]);
        }
    }

    int k = 0;
    for  (int i = 0;  i < N;  i ++) {
      if  (! A [i] . Good) {
        if  (i != k)
          A [k] = A [i];
        k ++;
      }
    }
    N = k;
  }  while  (N > 0);

  return count;
}


void ClusterMatches::Print_Cluster(const std::vector<std::vector<Match_t>>& clusters, const char* label) {
  if(clusters.empty()) {
    std::cout << label << '\n';
    return;
  }

  for(const auto cl : clusters) {
    std::cout << label << '\n'
              << std::setw(8) << cl[0].Start1 << ' '
              << std::setw(8) << cl[0].Start2 << ' '
              << std::setw(6) << cl[0].Len << ' '
              << "   none      -      -\n";
    for(size_t i = 1; i < cl.size(); ++i) {
      const long int adj = cl[i] . Simple_Adj;
      std::cout << std::setw(8) << (cl[i].Start1 + adj) << ' '
                << std::setw(8) << (cl[i].Start2 + adj) << ' '
                << std::setw(6) << (cl[i].Len - adj) << ' ';
      if  (adj == 0)
        std::cout << std::setw(7) << "none" << ' ';
      else
        std::cout << std::setw(7) << (-adj) << ' ';
      std::cout << std::setw(6) << (cl[i].Start1 + adj - cl[i-1].Start1 - cl[i-1].Len) << ' '
                << std::setw(6) << (cl[i].Start2 + adj - cl[i-1].Start2 - cl[i-1].Len) << '\n';
    }
    label = "#";
  }
}

int ClusterMatches::Process_Matches(Match_t * A, UnionFind& UF, int N, std::vector<std::vector<Match_t>>& clusters) {
  //  Process matches  A [1 .. N]  and output them after
  //  a line containing  label .

  //  Use Union-Find to create connected-components based on
  //  separation and similar diagonals between matches
  UF.reset(N);

  std::sort(A + 1, A + N + 1, By_Start2);
  Filter_Matches (A + 1, N);

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
    print_ct += Process_Cluster (A + i, cluster_size, clusters);
  }
  return print_ct;
}
} // namespace mgaps
} // namespace mummer
