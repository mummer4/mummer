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


#include "sw_alignscore.hh"
#include "tigrinc.hh"
#include <vector>
using namespace std;





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
struct Score
{
  long int value;
  char used;
};

struct Node
{
  Score S[3];
  Score * max;
};

struct Diagonal
{
  long int lbound, rbound;   // left(lower) and right(upper) bounds
  Node * I;          // the matrix nodes
};




//--------------------------------------------------------------- Externs ----//
extern int _break_len;
extern int _banding;
extern int _matrix_type;





//----------------------------------------- Private Function Declarations ----//
bool _alignEngine
     (const char * A0, long int Astart, long int & Aend,
      const char * B0, long int Bstart, long int & Bend,
      vector<long int> & Delta, unsigned int m_o);





//------------------------------------------- Public Function Definitions ----//
inline bool alignSearch
     (const char * A0, long int Astart, long int & Aend,
      const char * B0, long int Bstart, long int & Bend,
      unsigned int m_o)

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

{
  bool rv;
  vector<long int> n_v;

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

  rv = _alignEngine (A0, Astart, Aend, B0, Bstart, Bend, n_v, m_o);

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




inline bool alignTarget
     (const char * A0, long int Astart, long int & Aend,
      const char * B0, long int Bstart, long int & Bend,
      vector<long int> & Delta, unsigned int m_o)

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

  rv = _alignEngine (A0, Astart, Aend, B0, Bstart, Bend, Delta, m_o);

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




inline int getBreakLen
     ( )

     //  Returns the current value of _break_len

{
  return _break_len;
}

inline int getBanding()
{ return _banding; }



inline int getMatrixType
     ( )

     //  Returns the current value of _matrix_type

{
  return _matrix_type;
}




inline void setBreakLen
     (const int Break_Len)

     //  Resets the break length (see comments on DEFAULT_BREAK_LEN above)

{
  if ( Break_Len < 1  ||  Break_Len > MAX_ALIGNMENT_LENGTH )
    fprintf (stderr, "WARNING: Invalid break length %d, ignoring\n", Break_Len);
  else
    _break_len = Break_Len;
  return;
}

inline void setBanding(const int Banding)
{
  if ( Banding < 0 )
    fprintf (stderr, "WARNING: Invalid banding %d, ignoring\n", Banding);
  else
    _banding = Banding;
  return;
}



inline void setMatrixType
     (const int Matrix_Type)

     //  Resets the _matrix_type

{
  if ( Matrix_Type < 0 || Matrix_Type > 3 )
    fprintf (stderr,
	     "WARNING: Invalid matrix type %d, ignoring\n", Matrix_Type);
  else
    _matrix_type = Matrix_Type;
  return;
}



#endif // #ifndef __SW_ALIGN_HH
