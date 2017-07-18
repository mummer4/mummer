//------------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: sw_align.hh
//         Date: 7 / 18 / 2002
//
//   Description: Routines for DNA sequence alignment. See individual public
//               functions for a description of their applications. Designed
//              to be used by postNUC and postPRO for the extension of
//             alignments off of clusters identified by MUMmer.
//
//------------------------------------------------------------------------------

#ifndef __SW_ALIGN_HH
#define __SW_ALIGN_HH

//-- NOTE: these debug options will significantly hamper program performance
//#define _DEBUG_ASSERT      // performs assert functions to check validity
//#define _DEBUG_VERBOSE     // outputs various alignment statistics and values

#include <stdexcept>
#include <vector>

#include "sw_alignscore.hh"
#include "tigrinc.hh"

namespace mummer {
namespace sw_align {

//------------------------------------------------------------- Constants ----//
//-- Modus operandi bit masks 00001, 00010, 00100, 01000, 10000

static const unsigned int DIRECTION_BIT = 0x1;
static const unsigned int SEARCH_BIT = 0x2;
static const unsigned int FORCED_BIT = 0x4;
static const unsigned int OPTIMAL_BIT = 0x8;
static const unsigned int SEQEND_BIT = 0x10;


//-- Modus operandi of the primary alignment functions

static const unsigned int FORWARD_ALIGN = 0x1;
static const unsigned int OPTIMAL_FORWARD_ALIGN = 0x9;
static const unsigned int FORCED_FORWARD_ALIGN = 0x5;
//-- Align forward until target is reached OR score worsens
//   Fully uses memory and DOES create delta information
//   If "OPTIMAL" maximize the alignment score by shrinking coverage
//   If "FORCED" ignore score and force alignment to reach its target

static const unsigned int FORWARD_SEARCH = 0x3;
static const unsigned int OPTIMAL_FORWARD_SEARCH = 0xB;
static const unsigned int FORCED_FORWARD_SEARCH = 0x7;
//-- Align forward until target is reached OR score worsens
//   Only uses partial memory and DOES NOT create delta information
//   If "OPTIMAL" maximize the alignment score by shrinking coverage
//   If "FORCED" ignore score and force alignment to reach its target

static const unsigned int BACKWARD_SEARCH = 0x2;
static const unsigned int OPTIMAL_BACKWARD_SEARCH = 0xA;
static const unsigned int FORCED_BACKWARD_SEARCH = 0x6;
//-- Align backward until target is reached OR score worsens
//   Only uses partial memory and DOES NOT create delta information
//   If "OPTIMAL" maximize the alignment score by shrinking coverage
//   If "FORCED" ignore score and force alignment to reach its target

//-- Maximum number of bases (in either sequence) that the alignSearch may go
static const long int MAX_SEARCH_LENGTH = 10000;

//-- Maximum number of bases (in either sequence) that the alignTarget may go
static const long int MAX_ALIGNMENT_LENGTH = 10000;

//------------------------------------------------------ Type Definitions ----//
// An auto expanding vector
template<typename T>
class auto_vector {
  std::vector<T> m_vec;
public:
  typedef typename std::vector<T>               vector_type;
  typedef typename vector_type::value_type      value_type;
  typedef typename vector_type::reference       reference;
  typedef typename vector_type::const_reference const_reference;
  typedef typename vector_type::size_type       size_type;
  typedef typename vector_type::iterator        iterator;
  typedef typename vector_type::const_iterator  const_iterator;

  auto_vector() = default;

  reference operator[](size_type n) {
    if(__builtin_expect(n >= m_vec.size(), 0))
      m_vec.resize(n + 1);
    return m_vec[n];
  }
  const_reference operator[](size_type n) const {
    return m_vec[n];
  }
  void clear() noexcept {
    m_vec.clear();
  }
  void resize(size_type n) {
    m_vec.resize(n);
  }
  void resize(size_type n, const value_type& val) {
    m_vec.resize(n, val);
  }

  iterator begin() { return m_vec.begin(); }
  const_iterator begin() const { return m_vec.begin(); }
  const_iterator cbegin() const { return m_vec.cbegin(); }
  iterator end() { return m_vec.end(); }
  const_iterator end() const { return m_vec.end(); }
  const_iterator cend() const { return m_vec.cend(); }
};

struct Node
{
  long int values[3];
  char used[3];
  int m_max;

  long int max() const { return values[m_max]; }
  void max(int i) { m_max = i; }
  int edit() const { return m_max; }
};

struct Diagonal
{
  long int          lbound, rbound; // left(lower) and right(upper) bounds
  auto_vector<Node> I;          // the matrix nodes
};

// Auto expanding non-square matrix which minimizes allocation / free
class DiagonalMatrix {
  std::vector<Diagonal> m_diag;
  size_t                m_size; // Actual length.

public:
  DiagonalMatrix() : m_size(0) { }

  Diagonal& operator[](size_t n) {
    if(n >= m_size) {
      if(n >= m_diag.size())
        m_diag.resize(n + 1);
      m_size = n + 1;
    }
    return m_diag[n];
  }
  const Diagonal& operator[](size_t n) const {
    return m_diag[n];
  }
  void clear() noexcept {
    for(size_t i = 0; i < m_size; ++i)
      m_diag[i].I.clear();
    m_size = 0;
  }
};




class aligner {
  const int _break_len;
  const int _banding;
  const int _matrix_type;

public:
  aligner()
    : _break_len(200) // Number of bases to extend past global high score before giving up
    , _banding(0) // No banding by default
    , _matrix_type(NUCLEOTIDE)
  { }

  aligner(int break_len, int banding, int matrix_type)
    : _break_len(break_len)
    , _banding(banding)
    , _matrix_type(matrix_type)
  {
    if(break_len < 1 || break_len > MAX_ALIGNMENT_LENGTH)
      throw std::invalid_argument("Break length must be between 1 and MAX_ALIGNMENT_LENGTH included");
    if(banding < 0)
      throw std::invalid_argument("Banding must be >= 0");
    if(matrix_type < 0 || matrix_type > 3)
      throw std::invalid_argument("Matrix type must be between 0 and 3 included");
  }

  //------------------------------------------- Public Function Definitions ----//
  //  PURPOSE: This function aligns the sequences A0 and B0, starting
  //      at positions Astart and Bstart, as far as possible until the
  //      cumulative score does not improve for 'break_len' bases, or Aend
  //      and Bend are reached.  Aend and Bend are changed to reflect the
  //      ending position of the alignment, or if reached they stay the same.
  //      This function destroys the edit matrix as it goes along to save
  //      memory, therefore it only calculates the positions at which the
  //      alignment score falls off, it does not generate the alignment data
  //      (like as in alignTarget).
  //    INPUT: A0 and B0 are sequences such that A [1...N] and B [1...N].
  //      Usually the '\0' character is placed at the zero and end index
  //      of the arrays A and B. Astart and Bstart are the (inclusive)
  //      starting positions of the alignment. Aend and Bend are the
  //      (inclusive) target positions for the alignment. If the A / Bend
  //      positions are not reached, they are changed to the positions where
  //      the search terminated. These values are relative to the start of
  //      the sequences A and B. m_o is the modus operandi of the function,
  //      the m_o must be a Search, see beginning of file for descriptions.
  //   RETURN: This function returns 'true' if the alignment reached the
  //      target positions Aend and Bend, the function returns false if
  //      otherwise. Aend and Bend are changed to reflect the termination
  //      positions of the alignment.
  //  CONDITIONS:  A0 and B0 must not be null, and have a sequence length
  //      of atleast 1. Astart, Bstart, Aend and Bend must all be greater
  //      than zero and less or equal to the lengths of A and B (depending
  //      on which sequence they reference). If it is a forward extension
  //      Aend must be greater than Astart, and if it is a backward extension
  //      Astart must be greater than Aend. Astart must not equal Aend. The
  //      same rules apply for Bstart and Bend.
  inline bool alignSearch(const char * A0, long int Astart, long int & Aend,
                          const char * B0, long int Bstart, long int & Bend,
                          unsigned int m_o, DiagonalMatrix& Diag) const;
  inline bool alignSearch(const char * A0, long int Astart, long int & Aend,
                          const char * B0, long int Bstart, long int & Bend,
                          unsigned int m_o) const {
    DiagonalMatrix Diag;
    return alignSearch(A0, Astart, Aend,
                       B0, Bstart, Bend,
                       m_o, Diag);
  }

  //  PURPOSE: This function aligns the sequences A0 and B0, starting
  //      at positions Astart and Bstart, as far as possible until the
  //      cumulative score does not improve for 'break_len' bases, or Aend
  //      and Bend are reached.  Aend and Bend are changed to reflect the
  //      ending position of the alignment, or if reached they stay the same.
  //      This algorithm generates alignment data and does not destroy
  //      the edit matrix (like as in alignSearch).
  //    INPUT: A0 and B0 are sequences such that A [1...N] and B [1...N].
  //      Usually the '\0' character is placed at the zero and end index
  //      of the arrays A and B. Astart and Bstart are the (inclusive)
  //      starting positions of the alignment. Aend and Bend are the
  //      (inclusive) target positions for the alignment. If the A / Bend
  //      positions are not reached, they are changed to the positions where
  //      the search terminated. These values are relative to the start of
  //      the sequences A and B. Delta is an integer vector. It does not
  //      need to be empty, it will append the results to the end.
  //      m_o is the modus operandi of the function, the m_o must not be
  //      a search, and must be in the forward direction. See beginning of
  //      file for descriptions.
  //   RETURN: This function returns 'true' if the alignment reached the
  //      target positions Aend and Bend, the function returns false if
  //      otherwise. Aend and Bend are changed to reflect the termination
  //      positions of the alignment. The delta encoded alignment data
  //      is appended to the end of the Delta vector and terminated with
  //      a zero integer.
  //  CONDITIONS:  A0 and B0 must not be null, and have a sequence length
  //      of atleast 1. Astart, Bstart, Aend and Bend must all be greater
  //      than zero and less or equal to the lengths of A and B (depending
  //      on which sequence they reference). Aend must be greater than
  //      Astart. Astart must not equal Aend. The same rules apply for
  //      Bstart and Bend.
  inline bool alignTarget(const char * A0, long int Astart, long int & Aend,
                          const char * B0, long int Bstart, long int & Bend,
                          std::vector<long int>& Delta, unsigned int m_o,
                          DiagonalMatrix& Diag) const;
  inline bool alignTarget(const char * A0, long int Astart, long int & Aend,
                          const char * B0, long int Bstart, long int & Bend,
                          std::vector<long int>& Delta, unsigned int m_o) const {
    DiagonalMatrix Diag;
    return alignTarget(A0, Astart, Aend,
                       B0, Bstart, Bend,
                       Delta, m_o, Diag);
  }


  int breakLen() const { return _break_len; }
  int banding() const { return _banding; }
  int matrixType() const { return _matrix_type; }

  int good_score() const { return GOOD_SCORE[_matrix_type]; }
  int cont_gap_score() const { return CONT_GAP_SCORE[_matrix_type]; }
  int match_score(int i, int j) const { return MATCH_SCORE[_matrix_type][i][j]; }

protected:
  //----------------------------------------- Private Function Declarations ----//
  bool _alignEngine(const char * A0, long int Astart, long int & Aend,
                    const char * B0, long int Bstart, long int & Bend,
                    std::vector<long int> & Delta, unsigned int m_o, DiagonalMatrix& Diag) const;

  long int scoreMatch (const Diagonal& Diag, long int Dct, long int CDi,
                       const char * A, const char * B, long int N, unsigned int m_o) const;

};

// Identical to the above aligner class, with one difference. It keeps
// a DiagonalMatrix buffer around. This is a speed optimization
// (avoids repeated memory allocation/deallocation), but the
// alignSearch and alignTarget methods are not const anymore and are
// not thread safe.
class aligner_buffer : public aligner {
  mutable DiagonalMatrix m_Diag;
public:
  aligner_buffer() = default;
  aligner_buffer(int break_len, int banding, int matrix_type) : aligner(break_len, banding, matrix_type) { }

  // Warning: not thread safe!
  bool alignTarget(const char * A0, long int Astart, long int & Aend,
                   const char * B0, long int Bstart, long int & Bend,
                   std::vector<long int>& Delta, unsigned int m_o) const {
    return aligner::alignTarget(A0, Astart, Aend,
                                B0, Bstart, Bend,
                                Delta, m_o, m_Diag);
  }

  // Warning: not thread safe!
  bool alignSearch(const char * A0, long int Astart, long int & Aend,
                   const char * B0, long int Bstart, long int & Bend,
                   unsigned int m_o) const {
    return aligner::alignSearch(A0, Astart, Aend,
                                B0, Bstart, Bend,
                                m_o, m_Diag);
  }
};



bool aligner::alignSearch(const char * A0, long int Astart, long int & Aend,
                          const char * B0, long int Bstart, long int & Bend,
                          unsigned int m_o, DiagonalMatrix& Diag) const
{
  bool                  rv;
  std::vector<long int> n_v;

#ifdef _DEBUG_VERBOSE
  fprintf(stderr,"Running alignment search...\n");
  fprintf(stderr,"Astart = %ld, Atarget = %ld\n", Astart, Aend);
  fprintf(stderr,"Bstart = %ld, Btarget = %ld\n", Bstart, Bend);
  fprintf(stderr,"--------------------------------------\n");
#endif
#ifdef _DEBUG_ASSERT
  //-- Function pre-conditions
  assert ( A0 != NULL  &&  B0 != NULL );
  assert ( m_o & SEARCH_BIT );
  long int Alen = strlen ( A0 + 1 );
  long int Blen = strlen ( B0 + 1 );
  assert ( Astart > 0  &&  Aend > 0  &&  Astart <= Alen  &&  Aend <= Alen );
  assert ( Bstart > 0  &&  Bend > 0  &&  Bstart <= Blen  &&  Bend <= Blen );
  if ( m_o & DIRECTION_BIT )
    {
      assert ( Astart <= Aend  &&  Bstart <= Bend );
      assert ( Aend - Astart + 1 <= MAX_SEARCH_LENGTH  &&
	       Bend - Bstart + 1 <= MAX_SEARCH_LENGTH );
    }
  else
    {
      assert ( Astart >= Aend  &&  Bstart >= Bend );
      assert ( Astart - Aend + 1 <= MAX_SEARCH_LENGTH  &&
	       Bstart - Bend + 1 <= MAX_SEARCH_LENGTH);
    }
#endif

  rv = _alignEngine (A0, Astart, Aend, B0, Bstart, Bend, n_v, m_o, Diag);

#ifdef _DEBUG_VERBOSE
  fprintf(stderr,"--------------------------------------\n");
  fprintf(stderr,"Astart = %ld, Aend = %ld\n", Astart, Aend);
  fprintf(stderr,"Bstart = %ld, Bend = %ld\n", Bstart, Bend);
  fprintf(stderr,"Returned from alignment search\n");
#endif
#ifdef _DEBUG_ASSERT
  //-- Function post-conditions
  assert ( Aend > 0  &&  Aend <= Alen );
  assert ( Bend > 0  &&  Bend <= Blen );
  assert ( ~m_o & FORCED_BIT || rv == true );
#endif

  return rv;
}


bool aligner::alignTarget(const char * A0, long int Astart, long int & Aend,
                          const char * B0, long int Bstart, long int & Bend,
                          std::vector<long int> & Delta, unsigned int m_o,
                          DiagonalMatrix& Diag) const
{
  bool rv;

#ifdef _DEBUG_VERBOSE
  fprintf(stderr,"Running targeted alignment...\n");
  fprintf(stderr,"Astart = %ld, Atarget = %ld\n", Astart, Aend);
  fprintf(stderr,"Bstart = %ld, Btarget = %ld\n", Bstart, Bend);
  fprintf(stderr,"--------------------------------------\n");
#endif
#ifdef _DEBUG_ASSERT
  //-- Function pre-conditions
  assert ( m_o & DIRECTION_BIT  &&  ~m_o & SEARCH_BIT );
  assert ( A0 != NULL  &&  B0 != NULL );
  long int Alen = strlen ( A0 + 1 );
  long int Blen = strlen ( B0 + 1 );
  assert ( Astart > 0  &&  Aend > 0  &&  Astart <= Alen  &&  Aend <= Alen );
  assert ( Bstart > 0  &&  Bend > 0  &&  Bstart <= Blen  &&  Bend <= Blen );
  assert ( Astart <= Aend  &&  Bstart <= Bend );
  assert ( Aend - Astart + 1 <= MAX_ALIGNMENT_LENGTH  &&
	   Bend - Bstart + 1 <= MAX_ALIGNMENT_LENGTH);
#endif

  rv = _alignEngine (A0, Astart, Aend, B0, Bstart, Bend, Delta, m_o, Diag);

#ifdef _DEBUG_VERBOSE
  fprintf(stderr,"--------------------------------------\n");
  fprintf(stderr,"Astart = %ld, Aend = %ld\n", Astart, Aend);
  fprintf(stderr,"Bstart = %ld, Bend = %ld\n", Bstart, Bend);
  fprintf(stderr,"Returned from targeted alignment\n");
#endif
#ifdef _DEBUG_ASSERT
  //-- Function post-conditions
  assert ( Aend > 0  &&  Aend <= Alen );
  assert ( Bend > 0  &&  Bend <= Blen );
  assert ( ~m_o & FORCED_BIT || rv == true );
#endif

  return rv;
}

} // namespace sw_align
} // namespace mummer


#endif // #ifndef __SW_ALIGN_HH
