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

#include <mummer/tigrinc.hh>
#include <mummer/mgaps.hh>

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

int ClusterMatches::Filter_Matches(Match_t * A, const int N) {
//  Remove from  A [0 .. (N - 1)]  any matches that are internal to a repeat,
//  e.g., if seq1 has 27 As and seq2 has 20 then the first and
//  last matches will be kept, but the 6 matches in the middle will
//  be eliminated.  Also combine overlapping matches on the same
//  diagonal.  Pack all remaining matches into the front of  A  and
//  reduce the value of  N  if any matches are removed.
//  Matches in  A  *MUST* be sorted by  Start2  value.
//#pragma omp parallel for schedule(dynamic)
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

  const auto new_end = std::remove_if(A, A + N, [](const Match_t& m) { return !m.Good; });
  return new_end - A;
}


void ClusterMatches::Print_Cluster(const cluster_type& cl, const char* label, std::ostream& os) {
  os << label << '\n'
     << std::setw(8) << cl[0].Start1 << ' '
     << std::setw(8) << cl[0].Start2 << ' '
     << std::setw(6) << cl[0].Len << ' '
     << "   none      -      -\n";
  for(size_t i = 1; i < cl.size(); ++i) {
    const long int adj = cl[i] . Simple_Adj;
    os << std::setw(8) << (cl[i].Start1 + adj) << ' '
       << std::setw(8) << (cl[i].Start2 + adj) << ' '
       << std::setw(6) << (cl[i].Len - adj) << ' ';
    if  (adj == 0)
      os << std::setw(7) << "none" << ' ';
    else
      os << std::setw(7) << (-adj) << ' ';
    os << std::setw(6) << (cl[i].Start1 + adj - cl[i-1].Start1 - cl[i-1].Len) << ' '
       << std::setw(6) << (cl[i].Start2 + adj - cl[i-1].Start2 - cl[i-1].Len) << '\n';
  }
}

void ClusterMatches::Print_Clusters(const clusters_type& clusters, const char* label, std::ostream& os) {
  if(clusters.empty()) {
    os << label << '\n';
    return;
  }

  for(const auto& cl : clusters) {
    Print_Cluster(cl, label);
    label = "#";
  }
}

} // namespace mgaps
} // namespace mummer
