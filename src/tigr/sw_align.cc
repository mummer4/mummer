#include <math.h>
#include <limits>
#include <mummer/sw_align.hh>

namespace mummer {
namespace sw_align {

//-- Characters used in creating the alignment edit matrix, DO NOT ALTER!
static const int DELETE = 0;
static const int INSERT = 1;
static const int MATCH  = 2;
static const int START  = 3;
static const int NONE   = 4;

//----------------------------------------- Private Function Declarations ----//
static void generateDelta
     (const DiagonalMatrix& Diag, long int FinishCt, long int FinishCDi,
      long int N, std::vector<long int> & Delta);


static inline int maxScore(const Node & node);

static inline std::pair<long int, char> scoreEdit
     (const long int del, const long int ins, const long int mat);


//------------------------------------------ Private Function Definitions ----//
bool aligner::_alignEngine
     (const char* A0, long int Astart, long int & Aend,
      const char* B0, long int Bstart, long int & Bend,
      std::vector<long int>& Delta, unsigned int m_o,
      DiagonalMatrix& Diag) const

     //  A0 is a sequence such that A [1...\0]
     //  B0 is a sequence such that B [1...\0]
     //  The alignment should use bases A [Astart...Aend] (inclusive)
     //  The alignment should use bases B [Bstart...Bend] (inclusive)
     //       or [Aend...Astart] etc. if BACKWARD_SEARCH
     //       Aend must never equal Astart, same goes for Bend and Bstart
     //  Delta is an integer vector, not necessarily empty
     //  m_o is the modus operandi of the function:
     //      FORWARD_ALIGN, FORWARD_SEARCH, BACKWARD_SEARCH
     //  Returns true on success (Aend & Bend reached) or false on failure

{
  bool        TargetReached;    // the target was reached
  const char *A, *B;            // the sequence pointers to be used by this func

  static const long int min_score   = std::numeric_limits<long>::min(); // minimum possible score
  long int              high_score  = min_score; // global maximum score
  long int              xhigh_score = min_score; // non-optimal high score
  const long int        max_diff    = good_score() * _break_len; // max score difference

  long int CDi;                 // conceptual diagonal index (not relating to mem)
  long int Dct, Di;             // diagonal counter, actual diagonal index
  long int PDi, PPDi;           // previous diagonal index and prev prev diag index
  long int Ds, PDs, PPDs;       // diagonal size, prev, prev prev diagonal size where 'size' = rbound - lbound + 1

  long int Dl         = 2;      // current conceptual diagonal length
  long int lbound     = 0;      // current diagonal left(lower) node bound index
  long int rbound     = 0;      // current diagonal right(upper) node bound index
  long int FinishCt   = 0;      // diagonal containing the high_score
  long int FinishCDi  = 0;      // conceptual index of the high_score on FinishCt
  long int xFinishCt  = 0;      // non-optimal ...
  long int xFinishCDi = 0;      // non-optimal ...
  long int N, M, L;             // maximum matrix dimensions... N rows, M columns

  long int tlb, trb;
  double   Dmid  = .5;          // diag midpoint
  double   Dband = _banding/2.0; // diag banding

  int Iadj, Dadj, Madj;         // insert, delete and match adjust values

#ifdef _DEBUG_VERBOSE
  long int MaxL = 0;             // biggest diagonal seen
  long int TrimCt = 0;           // counter of nodes trimmed
  long int CalcCt = 0;           // counter of nodes calculated
#endif

  //-- Set up character pointers for the appropriate m_o
  if(m_o & DIRECTION_BIT) {
    A = A0 + ( Astart - 1 );
    B = B0 + ( Bstart - 1 );
    N = Aend - Astart + 1;
    M = Bend - Bstart + 1;
  } else {
    A = A0 + ( Astart + 1 );
    B = B0 + ( Bstart + 1 );
    N = Astart - Aend + 1;
    M = Bstart - Bend + 1;
  }

  //-- Initialize position 0,0 in the matrices
  Diag.clear();   // clear the list of diagonals to make up edit matrix
  { auto& D0 = Diag[0];
    D0.lbound = lbound;
    D0.rbound = rbound ++;
    auto& I0 = D0.I[0];
    I0.values[DELETE] = min_score;
    I0.values[INSERT] = min_score;
    I0.values[MATCH]  = 0;
    I0.max(MATCH);

    I0.used[DELETE] = NONE;
    I0.used[INSERT] = NONE;
    I0.used[MATCH]  = START;
  }

  L = N < M ? N : M;

  //-- **START** of diagonal processing loop
  //-- Calculate the rest of the diagonals until goal reached or score worsens
  for(Dct = 1; Dct <= N + M  && (Dct - FinishCt) <= _break_len  && lbound <= rbound; Dct++) {
    auto& CurD = Diag[Dct];
    CurD.lbound = lbound;
    CurD.rbound = rbound;

    //-- malloc space for the edit char and score nodes
    Ds = rbound - lbound + 1;

#ifdef _DEBUG_VERBOSE
    //-- Keep count of trimmed and calculated nodes
    CalcCt += Ds;
    TrimCt += Dl - Ds;
    if(Ds > MaxL )
      MaxL = Ds;
#endif

    //-- Set diagonal index adjustment values
    if(Dct <= N) {
      Iadj = 0;
      Madj = -1;
    } else {
      Iadj = 1;
      Madj = Dct == N + 1 ? 0 : 1;
    }
    Dadj = Iadj - 1;

    //-- Set parent diagonal values
    const auto& PrevD = Diag[Dct - 1] ; // previous diagonal
    PDs               = PrevD.rbound - PrevD.lbound + 1;
    PDi               = lbound + Dadj;
    PDi               = PDi - PrevD.lbound;

    //-- Set grandparent diagonal values
    const long PPDct = Dct - 2; //  prev prev diagonal
    const Diagonal& PPrevD = Diag[std::max((long)0, PPDct)]; // if PPDct < 0, not PPrevD is not used
    if(PPDct >= 0) {
      PPDs   = PPrevD.rbound - PPrevD.lbound + 1;
      PPDi   = lbound + Madj;
      PPDi   = PPDi - PPrevD.lbound;
    } else
      PPDi = PPDs = 0;

    //-- If forced alignment, don't keep track of global max
    if(m_o & FORCED_BIT )
      high_score = min_score;

    //-- **START** of internal node scoring loop
    //-- Calculate scores for every node (within bounds) for diagonal Dct
    for(CDi = lbound; CDi <= rbound; CDi ++) {
      //-- Set the index (in memory) of current node and clear score
      Di = CDi - CurD.lbound;
      auto& CurDi = CurD.I[Di];

      //-- Calculate DELETE score
      if(PDi >= 0  &&  PDi < PDs ) {
        const auto& PrevDi = PrevD.I[PDi];
        auto dummy = scoreEdit(
                  PrevDi.values[DELETE] + (PrevDi.used[DELETE] == NONE ? 0 : CONT_GAP_SCORE [_matrix_type]),
                  PrevDi.values[INSERT] + (PrevDi.used[INSERT] == NONE ? 0 : OPEN_GAP_SCORE [_matrix_type]),
                  PrevDi.values[MATCH]  + (PrevDi.used[MATCH] == NONE  ? 0 : OPEN_GAP_SCORE [_matrix_type]));
        CurDi.values[DELETE] = dummy.first;
        CurDi.used[DELETE] = dummy.second;
      } else {
        CurDi.values[DELETE] = min_score;
        CurDi.used[DELETE]  = NONE;
      }

      PDi ++;

      //-- Calculate INSERT score
      if(PDi >= 0  &&  PDi < PDs ) {
        const auto& PrevDi = PrevD.I[PDi];
        auto dummy = scoreEdit(
                  PrevDi.values[DELETE] + (PrevDi.used[DELETE] == NONE ? 0 : OPEN_GAP_SCORE [_matrix_type]),
                  PrevDi.values[INSERT] + (PrevDi.used[INSERT] == NONE ? 0 : CONT_GAP_SCORE [_matrix_type]),
                  PrevDi.values[MATCH]  + (PrevDi.used[MATCH] == NONE  ? 0 : OPEN_GAP_SCORE [_matrix_type]));
        CurDi.values[INSERT] = dummy.first;
        CurDi.used[INSERT] = dummy.second;
      } else {
        CurDi.values[INSERT] = min_score;
        CurDi.used[INSERT]  = NONE;
      }

      //-- Calculate MATCH/MIS-MATCH score
      if(PPDi >= 0  &&  PPDi < PPDs) {
        const auto& PPrevDi = PPrevD.I[PPDi];
        auto dummy = scoreEdit(
                  PPrevDi.values[DELETE],
                  PPrevDi.values[INSERT],
                  PPrevDi.values[MATCH]);
        CurDi.values[MATCH] = dummy.first;
        CurDi.used[MATCH] = dummy.second;
        CurDi.values[MATCH] += scoreMatch(CurD, Dct, CDi, A, B, N, m_o);
      } else {
        CurDi.values[MATCH] = min_score;
        CurDi.used[MATCH]  = NONE;
      }

      PPDi ++;

      CurDi.max(maxScore(CurDi));

      //-- Reset high_score if new global max was found
      if(CurDi.max() >= high_score) {
        high_score = CurDi.max();
        FinishCt   = Dct;
        FinishCDi  = CDi;
      }
    }
    //-- **END** of internal node scoring loop


    //-- Calculate max non-optimal score
    if(m_o & SEQEND_BIT  &&  Dct >= L) {
      if(L == N) {
        if(lbound == 0) {
          if(CurD.I[0].max() >= xhigh_score) {
            xhigh_score = CurD.I[0].max();
            xFinishCt   = Dct;
            xFinishCDi  = 0;
          }
        }
      } else  { // L == M
        if(rbound == M) {
          if(CurD.I[M-CurD.lbound].max() >= xhigh_score) {
            xhigh_score = CurD.I[M-CurD.lbound].max();
            xFinishCt   = Dct;
            xFinishCDi  = M;
          }
        }
      }
    }


    //-- If in extender modus operandi, free soon to be greatgrandparent diag
    // if(m_o & SEARCH_BIT  &&  Dct > 1 )
    //   free ( PPrevD->I );


    //-- Trim hopeless diagonal nodes
    const auto CurDi_start = CurD.I.cbegin();
    const auto CurDi_end   = CurDi_start + Ds;
    for(auto it = CurDi_start ; it < CurDi_end; ++it) {
      if(high_score - it->max() > max_diff )
        lbound ++;
      else
        break;
    }
    for(auto it = CurDi_end - 1; it >= CurDi_start; --it) {
      if(high_score - it->max() > max_diff )
        rbound --;
      else
        break;
    }

    //-- Grow new diagonal and reset boundaries
    if(Dct < N && Dct < M) {
      Dl ++; rbound ++; Dmid = (Dct+1)/2.0;
    } else if(Dct >= N && Dct >= M) {
      Dl --; lbound --; Dmid = N - (Dct+1)/2.0;
    } else if(Dct >= N) {
      lbound --; Dmid = N - (Dct+1)/2.0;
    } else {
      rbound ++; Dmid = (Dct+1)/2.0;
    }

    //-- Trim at hard band
    if(Dband > 0) {
      tlb = (long int)ceil(Dmid - Dband);
      if(lbound < tlb )
        lbound = tlb;
      trb = (long int)floor(Dmid + Dband);
      if(rbound > trb )
        rbound = trb;
    }

    if(lbound < 0 )
      lbound = 0;
    if(rbound >= Dl )
      rbound = Dl - 1;
  }
  //-- **END** of diagonal processing loop
  Dct --;

  //-- Check if the target was reached
  //   If OPTIMAL, backtrack to last high_score to maximize alignment score
  TargetReached = false;
  if(Dct == N + M) {
    if(~m_o & OPTIMAL_BIT || m_o & SEQEND_BIT) {
      TargetReached = true;
      FinishCt      = N + M;
      FinishCDi     = 0;
    } else if(FinishCt == Dct )
      TargetReached = true;
  } else if(m_o & SEQEND_BIT  &&  xFinishCt != 0) {
    //-- non-optimal, extend alignment to end of shortest seq if possible
    FinishCt  = xFinishCt;
    FinishCDi = xFinishCDi;
  }

  //-- Set A/Bend to finish positions
  long int Aadj = FinishCt <= N ? FinishCt - FinishCDi - 1 : N - FinishCDi - 1;
  long int Badj = FinishCt <= N ? FinishCDi - 1 : FinishCt - N + FinishCDi - 1;
  if(~m_o & DIRECTION_BIT) {
    Aadj *= -1;
    Badj *= -1;
  }
  Aend = Astart + Aadj;
  Bend = Bstart + Badj;

#ifdef _DEBUG_VERBOSE
  assert (FinishCt > 1);

  //-- Ouput calculation statistics
  if(TargetReached )
    fprintf(stderr,"Finish score = %ld : %ld,%ld\n",
	    Diag[FinishCt].I[0].max(), N, M);
  else
    fprintf(stderr,"High score = %ld : %ld,%ld\n", high_score,
	    labs(Aadj) + 1, labs(Badj) + 1);
  fprintf(stderr, "%ld nodes calculated, %ld nodes trimmed\n", CalcCt, TrimCt);
  if(m_o & DIRECTION_BIT )
    fprintf(stderr, "%ld bytes used\n",
	    (long int)sizeof(Diagonal) * Dct + (long int)sizeof(Node) * CalcCt);
  else
    fprintf(stderr, "%ld bytes used\n",
	    ((long int)sizeof(Diagonal) + (long int)sizeof(Node) * MaxL) * 2);
#endif

  //-- If in forward alignment m_o, create the Delta information
  if(~m_o & SEARCH_BIT )
    generateDelta(Diag, FinishCt, FinishCDi, N, Delta);

  return TargetReached;
}

long int aligner::scoreMatch
     (const Diagonal& Diag, long int Dct, long int CDi,
      const char * A, const char * B, long int N, unsigned int m_o) const

     //  Diag is the single diagonal that contains the node to be scored
     //  Dct is Diag's diagonal index in the edit matrix
     //  CDi is the conceptual node to be scored in Diag
     //  A and B are the alignment sequences
     //  N is the alignment target index in A
     //  m_o is the modus operandi of the alignment:
     //      FORWARD_ALIGN, FORWARD_SEARCH, BACKWARD_SEARCH

{
  int Dir;
  char Ac, Bc;

  //-- 1 for forward, -1 for reverse
  Dir = m_o & DIRECTION_BIT ? 1 : -1;

  //-- Locate the characters that need to be compared
  if ( Dct <= N )
    {
      Ac = *( A + ( (Dct - CDi) * Dir ) );
      Bc = *( B + ( (CDi) * Dir ) );
    }
  else
    {
      Ac = *( A + ( (N - CDi) * Dir ) );
      Bc = *( B + ( (Dct - N + CDi) * Dir ) );
    }

  if ( ! isalpha(Ac) )
    Ac = STOP_CHAR;
  if ( ! isalpha(Bc) )
    Bc = STOP_CHAR;

  return MATCH_SCORE [_matrix_type] [toupper(Ac) - 'A'] [toupper(Bc) - 'A'];
}



static void generateDelta
     (const DiagonalMatrix& Diag, long int FinishCt, long int FinishCDi,
      long int N, std::vector<long int> & Delta)

     //  Diag is the list of diagonals that compose the edit matrix
     //  FinishCt is the diagonal that contains the finishing node
     //  FinishCDi is the conceptual finishing node, in FinishCt, for the align
     //  N & M are the target positions for the alignment
     //  Delta is the vector in which to store the alignment data, new data
     //      will be appended onto any existing data.
     //  NOTE: there will be no zero at the end of the data, end of data
     //        is signaled by the end of the vector
     //  Return is void

{
  //-- Function pre-conditions
#ifdef _DEBUG_ASSERT
  assert ( Diag != NULL );
  assert ( FinishCt > 1 );
#endif

  long int Count;                // delta counter
  long int Dct = FinishCt;  // diagonal index
  long int CDi = FinishCDi; // conceptual node index
  long int Di = 0;          // actual node index
  long int Pi = 0;          // path index
  long int PSize = 100;     // capacity of the path space
  char * Reverse_Path;       // path space

  int edit;

  //-- malloc space for the edit path
  Reverse_Path = (char *) Safe_malloc ( PSize * sizeof(char) );

  //-- Which Score index is the maximum value in? Store in edit
  Di = CDi - Diag[Dct] . lbound;
  edit = Diag[Dct].I[Di].edit();


  //-- Walk the path backwards through the edit space
  while ( Dct >= 0 ) {
    //-- remalloc path space if neccessary
    if ( Pi >= PSize ) {
      PSize *= 2;
      Reverse_Path = (char *) Safe_realloc( Reverse_Path, sizeof(char) * PSize );
    }

    Di = CDi - Diag[Dct].lbound;
    auto used = Diag[Dct].I[Di].used[edit];

    Reverse_Path[Pi ++] = edit;
    switch ( edit ) {
    case DELETE :
      CDi = Dct -- <= N ? CDi - 1 : CDi;
      break;
    case INSERT :
      CDi = Dct -- <= N ? CDi : CDi + 1;
      break;
    case MATCH :
      CDi = Dct <= N ? CDi - 1 : ( Dct == N + 1 ? CDi : CDi + 1 );
      Dct -= 2;
      break;
    case START :
      Dct = -1;
      break;
    default :
      fprintf(stderr,"\nERROR: Invalid edit matrix entry,\n"
              "       please file a bug report\n");
      exit ( EXIT_FAILURE );
    }

    edit = used;
  }

  //-- Generate the delta information
  Count = 1;
  for (Pi -= 2; Pi >= 0; Pi --) {
    switch ( Reverse_Path[Pi] ) {
    case DELETE :
      Delta.push_back(-Count);
      Count = 1;
      break;
    case INSERT :
      Delta.push_back(Count);
      Count = 1;
      break;
    case MATCH :
      Count ++;
      break;
    case START :
      break;
    default :
      fprintf(stderr,"\nERROR: Invalid path matrix entry,\n"
              "       please file a bug report\n");
      exit ( EXIT_FAILURE );
    }
  }

  free (Reverse_Path);

  return;
}




static inline int maxScore(const Node & node)

     //  Return a pointer to the maximum score in the score array

{
  const auto& values = node.values;
  if(values[DELETE] > values[INSERT])
    return values[DELETE] > values[MATCH] ? DELETE : MATCH;
  else
    return values[INSERT] > values[MATCH] ? INSERT : MATCH;
}




static inline std::pair<long int, char> scoreEdit
     (const long int del, const long int ins, const long int mat)

     //  Assign current edit a maximal score using either del, ins or mat

{
  if ( del > ins ) {
    if ( del > mat ) {
      return std::make_pair(del, DELETE);
    } else {
      return std::make_pair(mat, MATCH);
    }
  } else if ( ins > mat ) {
    return std::make_pair(ins, INSERT);
  }  else {
    return std::make_pair(mat, MATCH);
  }
}

} // namespace sw_align
} // namespace mummer
