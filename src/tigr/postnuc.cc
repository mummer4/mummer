
//-- NOTE: this option will significantly hamper program performance,
//         mostly the alignment extension performance (sw_align.h)
//#define _DEBUG_ASSERT       // self testing assert functions

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>

#include <mummer/postnuc.hh>
#include <mummer/tigrinc.hh>
#include <mummer/sw_align.hh>

namespace mummer {
namespace postnuc {
std::ostream& operator<<(std::ostream& os, const Alignment& al) {
  os << '<' << al.sA << '-' << al.eA
     << " | " << al.sB << '-' << al.eB
     << " | " << al.Errors
     << " |";
  for(long x : al.delta)
    os << ' ' << x;
  return os << '>';
}

error_iterator_type& error_iterator_type::operator++()  {
  switch(m_error.type) {
  case NONE:
  case MISMATCH:
    ++m_error.dst;
    ++m_error.ref;
    m_error.qry += m_al.dirB;
    break;

  case INSERTION:
    m_error.dst = 1;
    ++m_k;
    ++m_error.ref;
    break;

  case DELETION:
    m_error.dst = 1;
    ++m_k;
    m_error.qry += m_al.dirB;
    break;
  }

  while(m_error.ref < m_ref_end) {
    if(m_k >= m_al.delta.size() || m_error.dst != std::abs(m_al.delta[m_k])) { // not at an indel
      if(*m_error.ref != (m_al.dirB == 1 ? *m_error.qry : comp(*m_error.qry))) {
        m_error.type = MISMATCH;
        break;
      }
      ++m_error.dst;
      ++m_error.ref;
      m_error.qry += m_al.dirB;
    } else { // indel
      m_error.type = m_al.delta[m_k] > 0 ? INSERTION : DELETION;
      break;
    }
  }
  return *this;
}

// Read one sequence from fasta stream into T (starting at index 1, 0
// is unused), store its name and returns true if successful.
bool Read_Sequence(std::istream& is, std::string& T, std::string& name) {
  int c = is.peek();
  for( ; c != EOF && c != '>'; c = is.peek())
    ignore_line(is);
  if(c == EOF) return false;
  std::getline(is, name);
  name = name.substr(1, name.find_first_of(" \t\n") - 1);

  T = '\0';
  for(c = is.peek(); c != EOF && c != '>'; c = is.peek()) {
    c = is.get();
    if(isspace(c)) continue;
    c = tolower(c);
    if(!isalpha(c) && c != '*')
      c = 'x';
    T += c;
  }
  return true;
}

void merge_syntenys::extendClusters(std::vector<Cluster> & Clusters,
                                    const char* Aseq, const long Alen, const char* Bseq, const long Blen,
                                    std::vector<Alignment>& Alignments /* the vector of alignment objects */) const

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

  const char* const A = Aseq;
  const char*       B;          // the sequences A and B
  std::unique_ptr<char[]> Brev; // the reverse complement of B

  unsigned int m_o;
  long int targetA, targetB;           // alignment extension targets in A and B

  std::vector<Match>::iterator Mp;          // match pointer

  std::vector<Alignment>::iterator CurrAp = Alignments.begin( );   // current align
  std::vector<Alignment>::iterator TargetAp;                // target align


  //-- Extend each cluster
  auto TargetCp = Clusters.end(); // the target cluster
  auto PrevCp   = Clusters.begin( ); // where the extensions last left off
  auto CurrCp   = Clusters.begin( ); // the current cluster being extended
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
      B = Bseq;
    else if ( Brev )
      B = Brev.get();
    else {
      Brev.reset(new char[Blen + 2]);
      memcpy ( Brev.get() + 1, Bseq + 1, Blen );
      Brev[0] = Brev[Blen+1] = '\0';
      Reverse_Complement (Brev.get(), 1, Blen);
      B = Brev.get();
    }

    //-- Extend each match in the cluster
    for ( Mp = CurrCp->matches.begin( ); Mp < CurrCp->matches.end( ); ++Mp) {
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
        assert(CurrAp->sA >= 1 && CurrAp->eA <= Alen);
        assert(CurrAp->sB >= 1 && CurrAp->eB <= Blen);
      } else { //-- Create a new alignment object
        Alignments.push_back({ *Mp, CurrCp->dirB } );
        CurrAp = Alignments.end( ) - 1;

        if ( DO_EXTEND  ||  Mp != CurrCp->matches.begin ( ) ) {
          //-- Target the closest/best alignment object
          TargetAp = getReverseTargetAlignment (Alignments, CurrAp);
          assert(CurrAp->sA >= 1 && CurrAp->eA <= Alen);
          assert(CurrAp->sB >= 1 && CurrAp->eB <= Blen);
          assert(TargetAp >= Alignments.begin());
          assert(TargetAp <= Alignments.end());

          //-- Extend the new alignment object backwards
          if ( extendBackward (Alignments, CurrAp, TargetAp, A, B) )
            CurrAp = TargetAp;
          assert(CurrAp->sA >= 1 && CurrAp->eA <= Alen);
          assert(CurrAp->sB >= 1 && CurrAp->eB <= Blen);
        }
      }

      m_o = sw_align::FORWARD_ALIGN;

      //-- Try to extend the current match forwards
      if ( Mp < CurrCp->matches.end( ) - 1 ) {
        //-- Target the next match in the cluster
        targetA = (Mp + 1)->sA;
        targetB = (Mp + 1)->sB;

        //-- Extend the current alignment object forward
        target_reached = extendForward (CurrAp, A, targetA, B, targetB, m_o);
      } else if ( DO_EXTEND ) {
        targetA = Alen;
        targetB = Blen;

        //-- Target the closest/best match in a future cluster
        TargetCp = getForwardTargetCluster (Clusters, CurrCp, targetA, targetB);
        assert(targetA <= Alen);
        assert(targetB <= Blen);
        if ( TargetCp == Clusters.end( ) ) {
          m_o |= sw_align::OPTIMAL_BIT;
          if ( TO_SEQEND )
            m_o |= sw_align::SEQEND_BIT;
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
  validateData (Alignments, Clusters, Aseq, Alen, Bseq, Blen);
#endif

  //-- Generate the error counts
  parseDelta(Alignments, Aseq, Bseq, Blen);
}

bool merge_syntenys::extendBackward(std::vector<Alignment> & Alignments, std::vector<Alignment>::iterator CurrAp,
                                    std::vector<Alignment>::iterator TargetAp, const char * A, const char * B) const

//  Extend an alignment backwards off of the current alignment object.
//  The current alignment object must be freshly created and consist
//  only of an exact match (i.e. the delta vector MUST be empty).
//  If the TargetAp alignment object is reached by the extension, it will
//  be merged with CurrAp and CurrAp will be destroyed. If TargetAp is
//  NULL the function will extend as far as possible. It is a strange
//  and dangerous function because it can delete CurrAp, so edit with
//  caution. Returns true if TargetAp was reached and merged, else false
//  Designed only as a subroutine for extendClusters, should be used
//  nowhere else.

{
  bool target_reached = false;
  bool overflow_flag  = false;
  bool double_flag    = false;

  std::vector<long int>::iterator Dp;

  unsigned int m_o;
  long int targetA, targetB;

  m_o = sw_align::BACKWARD_SEARCH;

  //-- Set the target coordinates
  if ( TargetAp != Alignments.end( ) )
    {
      targetA = TargetAp->eA;
      targetB = TargetAp->eB;
    }
  else
    {
      targetA = 1;
      targetB = 1;
      m_o |= sw_align::OPTIMAL_BIT;
    }

  //-- If alignment is too long, bring with bounds and set overflow_flag true
  if ( CurrAp->sA - targetA + 1 > sw_align::MAX_ALIGNMENT_LENGTH )
    {
      targetA = CurrAp->sA - sw_align::MAX_ALIGNMENT_LENGTH + 1;
      overflow_flag = true;
      m_o |= sw_align::OPTIMAL_BIT;
    }
  if ( CurrAp->sB - targetB + 1 > sw_align::MAX_ALIGNMENT_LENGTH )
    {
      targetB = CurrAp->sB - sw_align::MAX_ALIGNMENT_LENGTH + 1;
      if ( overflow_flag )
        double_flag = true;
      else
        overflow_flag = true;
      m_o |= sw_align::OPTIMAL_BIT;
    }


  if ( TO_SEQEND && !double_flag )
    m_o |= sw_align::SEQEND_BIT;


  //-- Attempt to reach the target
  target_reached = aligner.alignSearch (A, CurrAp->sA, targetA,
                                B, CurrAp->sB, targetB, m_o);

  if ( overflow_flag  ||  TargetAp == Alignments.end( ) )
    target_reached = false;

  if ( target_reached )
    {
      //-- Merge the two alignment objects
      extendForward (TargetAp, A, CurrAp->sA,
                     B, CurrAp->sB, sw_align::FORCED_FORWARD_ALIGN);
      TargetAp->eA = CurrAp->eA;
      TargetAp->eB = CurrAp->eB;
      Alignments.pop_back( );
    }
  else
    {
      aligner.alignTarget (A, targetA, CurrAp->sA,
                           B, targetB, CurrAp->sB,
                           CurrAp->delta, sw_align::FORCED_FORWARD_ALIGN);
      CurrAp->sA = targetA;
      CurrAp->sB = targetB;

      //-- Update the deltaApos value for the alignment object
      for ( Dp = CurrAp->delta.begin( ); Dp < CurrAp->delta.end( ); Dp ++ )
        CurrAp->deltaApos += *Dp > 0 ? *Dp : labs(*Dp) - 1;
    }

  return target_reached;
}




bool merge_syntenys::extendForward(std::vector<Alignment>::iterator CurrAp, const char * A, long int targetA,
                                   const char * B, long int targetB, unsigned int m_o) const

//  Extend an alignment forwards off the current alignment object until
//  target or end of sequence is reached, and merge the delta values of the
//  alignment object with the new delta values generated by the extension.
//  Return true if the target was reached, else false

{
  long int ValA;
  long int ValB;
  unsigned int Di;
  bool target_reached;
  bool overflow_flag = false;
  bool double_flag = false;
  std::vector<long int>::iterator Dp;

  //-- Set Di to the end of the delta vector
  Di = CurrAp->delta.size( );

  //-- Calculate the distance between the start and end positions
  ValA = targetA - CurrAp->eA + 1;
  ValB = targetB - CurrAp->eB + 1;

  //-- If the distance is too long, shrink it and set the overflow_flag
  if ( ValA > sw_align::MAX_ALIGNMENT_LENGTH )
    {
      targetA = CurrAp->eA + sw_align::MAX_ALIGNMENT_LENGTH - 1;
      overflow_flag = true;
      m_o |= sw_align::OPTIMAL_BIT;
    }
  if ( ValB > sw_align::MAX_ALIGNMENT_LENGTH )
    {
      targetB = CurrAp->eB + sw_align::MAX_ALIGNMENT_LENGTH - 1;
      if ( overflow_flag )
        double_flag = true;
      else
        overflow_flag = true;
      m_o |= sw_align::OPTIMAL_BIT;
    }

  if ( double_flag )
    m_o &= ~sw_align::SEQEND_BIT;

  target_reached = aligner.alignTarget (A, CurrAp->eA, targetA,
                                        B, CurrAp->eB, targetB,
                                        CurrAp->delta, m_o);

  //-- Notify user if alignment was chopped short
  if ( target_reached  &&  overflow_flag )
    target_reached = false;

  //-- Pick up delta where left off (Di) and merge with new deltas
  if ( Di < CurrAp->delta.size( ) )
    {
      //-- Merge the deltas
      ValA = (CurrAp->eA - CurrAp->sA + 1) - CurrAp->deltaApos - 1;
      CurrAp->delta[Di] += CurrAp->delta[Di] > 0 ? ValA : -(ValA);
      if ( CurrAp->delta[Di] == 0  ||  ValA < 0 ) {
        std::cerr << "ERROR: failed to merge alignments at position " << CurrAp->eA << '\n'
             << "       Please file a bug report\n";
        exit (EXIT_FAILURE);
      }

      //-- Update the deltaApos
      for ( Dp = CurrAp->delta.begin( ) + Di; Dp < CurrAp->delta.end( ); Dp ++ )
        CurrAp->deltaApos += *Dp > 0 ? *Dp : labs(*Dp) - 1;
    }

  //-- Set the alignment coordinates
  CurrAp->eA = targetA;
  CurrAp->eB = targetB;

  return target_reached;
}

std::vector<Cluster>::iterator merge_syntenys::getForwardTargetCluster
(std::vector<Cluster> & Clusters, std::vector<Cluster>::iterator CurrCp,
 long int & targetA, long int & targetB) const

//  Return the cluster that is most likely to successfully join (in a
//  forward direction) with the current cluster. The returned cluster
//  must contain 1 or more matches that are strictly greater than the end
//  of the current cluster. The targeted cluster must also be on a
//  diagonal close enough to the current cluster, so that a connection
//  could possibly be made by the alignment extender. Assumes clusters
//  have been sorted via AscendingClusterSort. Returns targeted cluster
//  and stores the target coordinates in targetA and targetB. If no
//  suitable cluster was found, the function will return NULL and target
//  A and targetB will remain unchanged.

{
  std::vector<Match>::iterator Mip;               // match iteratrive pointer
  std::vector<Cluster>::iterator Cp;              // cluster pointer
  std::vector<Cluster>::iterator Cip;             // cluster iterative pointer
  long int eA, eB;                           // possible target
  long int greater, lesser;                  // gap sizes between two clusters
  long int sA = CurrCp->matches.rbegin( )->sA +
    CurrCp->matches.rbegin( )->len - 1;      // the endA of the current cluster
  long int sB = CurrCp->matches.rbegin( )->sB +
    CurrCp->matches.rbegin( )->len - 1;      // the endB of the current cluster

  //-- End of sequences is the default target, set distance accordingly
  long int dist = (targetA - sA < targetB - sB ? targetA - sA : targetB - sB);

  //-- For all clusters greater than the current cluster (on sequence A)
  Cp = Clusters.end( );
  for ( Cip = CurrCp + 1; Cip < Clusters.end( ); Cip ++ )
    {
      //-- If the cluster is on the same direction
      if ( CurrCp->dirB == Cip->dirB )
        {
          eA = Cip->matches.begin( )->sA;
          eB = Cip->matches.begin( )->sB;

          //-- If the cluster overlaps the current cluster, strip some matches
          if ( ( eA < sA  ||  eB < sB )  &&
               Cip->matches.rbegin( )->sA >= sA  &&
               Cip->matches.rbegin( )->sB >= sB )
            {
              for ( Mip = Cip->matches.begin( );
                    Mip < Cip->matches.end( )  &&  ( eA < sA  ||  eB < sB );
                    Mip ++ )
                {
                  eA = Mip->sA;
                  eB = Mip->sB;
                }
            }

          //-- If the cluster is strictly greater than current cluster
          if ( eA >= sA  &&  eB >= sB )
            {
              if ( eA - sA > eB - sB )
                {
                  greater = eA - sA;
                  lesser = eB - sB;
                }
              else
                {
                  lesser = eA - sA;
                  greater = eB - sB;
                }

              //-- If the cluster is close enough
              if ( greater < aligner.breakLen( )  ||
                   (lesser) * aligner.good_score() +
                   (greater - lesser) * aligner.cont_gap_score() >= 0 )
                {
                  Cp = Cip;
                  targetA = eA;
                  targetB = eB;
                  break;
                }
              else if ( (greater << 1) - lesser < dist )
                {
                  Cp = Cip;
                  targetA = eA;
                  targetB = eB;
                  dist = (greater << 1) - lesser;
                }
            }
        }
    }


  return Cp;
}




std::vector<Alignment>::iterator merge_syntenys::getReverseTargetAlignment
(std::vector<Alignment> & Alignments, std::vector<Alignment>::iterator CurrAp) const

//  Return the alignment that is most likely to successfully join (in a
//  reverse direction) with the current alignment. The returned alignment
//  must be strictly less than the current cluster and be on a diagonal
//  close enough to the current alignment, so that a connection
//  could possibly be made by the alignment extender. Assumes clusters
//  have been sorted via AscendingClusterSort and processed in order, so
//  therefore all alignments are in order by their start A coordinate.

{
  long int greater, lesser;              // gap sizes between the two alignments
  const long int sA = CurrAp->sA;              // the startA of the current alignment
  const long int sB = CurrAp->sB;              // the startB of the current alignment

  //-- Beginning of sequences is the default target, set distance accordingly
  long int dist = (sA < sB ? sA : sB);

  //-- For all alignments less than the current alignment (on sequence A)
  auto Ap = Alignments.end( ); // alignment pointer
  for (auto Aip = CurrAp - 1; Aip >= Alignments.begin( ); --Aip) {
    //-- If the alignment is on the same direction
    if ( CurrAp->dirB == Aip->dirB ) {
      const long eA = Aip->eA;
      const long eB = Aip->eB;

      //-- If the alignment is strictly less than current cluster
      if ( eA <= sA  && eB <= sB ) {
        if ( sA - eA > sB - eB ) {
          greater = sA - eA;
          lesser  = sB - eB;
        } else {
          lesser  = sA - eA;
          greater = sB - eB;
        }

        //-- If the cluster is close enough
        if ( greater < aligner.breakLen( )  ||
             (lesser) * aligner.good_score() + (greater - lesser) * aligner.cont_gap_score() >= 0 ) {
          Ap = Aip;
          break;
        } else if ( (greater << 1) - lesser < dist ) {
          Ap = Aip;
          dist = (greater << 1) - lesser;
        }
      }
    }
  }

  return Ap;
}




bool isShadowedCluster(std::vector<Cluster>::const_iterator CurrCp,
                       const std::vector<Alignment> & Alignments, std::vector<Alignment>::const_iterator Ap)

//  Check if the current cluster is shadowed by a previously produced
//  alignment region. Return true if it is, else false.

{
  const long int sA = CurrCp->matches.front().sA;
  const long int eA = CurrCp->matches.back().sA + CurrCp->matches.back().len - 1;
  const long int sB = CurrCp->matches.front().sB;
  const long int eB = CurrCp->matches.back().sB + CurrCp->matches.back().len - 1;

  if ( ! Alignments.empty( ) ) {           // if there are alignments to use
    //-- Look backwards in hope of finding a shadowing alignment
    for(auto Aip = Ap ; Aip != Alignments.cbegin( ); -- Aip ) {
      //-- If in the same direction and shadowing the current cluster, break
      if ( Aip->dirB == CurrCp->dirB &&
           Aip->eA >= eA  &&  Aip->eB >= eB &&
           Aip->sA <= sA  &&  Aip->sB <= sB )
        return true; // shadow found
    }
  }

  //-- Return false if Alignments was empty or loop was not broken
  return false;
}




void __parseAbort
(const char * s, const char* file, size_t line)

//  Abort the program if there was an error in parsing file 's'

{
  std::cerr << "ERROR: " << file << ':' << line
            << " Could not parse input from '" << s << "'. \n"
            << "Please check the filename and format, or file a bug report\n";
  exit (EXIT_FAILURE);
}

void merge_syntenys::parseDelta
(std::vector<Alignment> & Alignments,
 const char* const Aseq, const char* const Bseq, const long Blen) const

// Use the delta information to generate the error counts for each
// alignment, and fill this information into the data type

{
  const char* const A = Aseq;
  const char* B;
  std::unique_ptr<char[]> Brev;
  char ch1, ch2;
  long int Delta;
  int Sign;
  long int i, Apos, Bpos;
  long int Remain, Total;
  long int Errors, SimErrors;
  long int NonAlphas;
  std::vector<Alignment>::iterator Ap;
  std::vector<long int>::iterator Dp;

  for ( Ap = Alignments.begin( ); Ap != Alignments.end( ); ++Ap) {
    B = Bseq;
      if ( Ap->dirB == REVERSE_CHAR ) {
        if (!Brev) {
          Brev.reset(new char[Blen + 2]);
          //          strcpy ( Brev.get() + 1, B + 1 );
          memcpy(Brev.get() + 1, B + 1, Blen);
          Brev[0] = Brev[Blen+1] = '\0';
          Reverse_Complement (Brev.get(), 1, Blen);
        }
        B = Brev.get();
      }

      Apos = Ap->sA;
      Bpos = Ap->sB;

      Errors = 0;
      SimErrors = 0;
      NonAlphas = 0;
      Remain = Ap->eA - Ap->sA + 1;
      Total = Remain;

      //-- For all delta's in this alignment
      for ( Dp = Ap->delta.begin( ); Dp != Ap->delta.end( ); ++Dp) {
        Delta = *Dp;
        Sign = Delta > 0 ? 1 : -1;
        Delta = std::abs ( Delta );

        //-- For all the bases before the next indel
          for ( i = 1; i < Delta; i ++ ) {
            ch1 = A [Apos ++];
            ch2 = B [Bpos ++];

            if ( !isalpha (ch1) ) {
              ch1 = sw_align::STOP_CHAR;
              NonAlphas ++;
            }
            if ( !isalpha (ch2) ) {
              ch2 = sw_align::STOP_CHAR;
              NonAlphas ++;
            }

            ch1 = toupper(ch1);
            ch2 = toupper(ch2);
            if (1 > aligner.match_score(ch1 - 'A', ch2 - 'A'))
              SimErrors ++;
            if ( ch1 != ch2 )
                Errors ++;
          }

          //-- Process the current indel
          Remain -= i - 1;
          Errors ++;
          SimErrors ++;

          if ( Sign == 1 ) {
            Apos ++;
            Remain --;
          } else {
            Bpos ++;
            Total ++;
          }
      }

      //-- For all the bases after the final indel
      for ( i = 0; i < Remain; i ++ ) {
        //-- Score character match and update error counters
        ch1 = A [Apos ++];
        ch2 = B [Bpos ++];

        if ( !isalpha (ch1) ) {
          ch1 = sw_align::STOP_CHAR;
          NonAlphas ++;
        }
        if ( !isalpha (ch2) ) {
          ch2 = sw_align::STOP_CHAR;
          NonAlphas ++;
        }

        ch1 = toupper(ch1);
        ch2 = toupper(ch2);
        if ( 1 > aligner.match_score(ch1 - 'A', ch2 - 'A') ) // 
          SimErrors ++;
        if ( ch1 != ch2 )
          Errors ++;
      }

      Ap->Errors = Errors;
      Ap->SimErrors = SimErrors;
      Ap->NonAlphas = NonAlphas;
  }
}


// always_assert: similar to assert macro, but not subject to NDEBUG
#define always_assert(x)                                                \
  if(!(x)) {                                                            \
    std::cerr << __FILE__ << ':' << __LINE__                            \
              << ": assertion failed " << #x << std::endl;              \
      abort();                                                          \
  }

void validateData
(std::vector<Alignment> Alignments, std::vector<Cluster> Clusters,
 const char* Aseq, const long Alen, const char* Bseq, const long Blen)

//  Self test function to check that the cluster and alignment information
//  is valid

{
  std::unique_ptr<char[]> Brev;
  std::vector<Cluster>::iterator   Cp;
  std::vector<Match>::iterator     Mp;
  std::vector<Alignment>::iterator Ap;
  const char* const                A = Aseq;
  const char*                      B;

  for ( Cp = Clusters.begin( ); Cp < Clusters.end( ); Cp ++ ) {
    always_assert ( Cp->wasFused );

    //-- Pick the right directional sequence for B
    if ( Cp->dirB == FORWARD_CHAR ) {
      B = Bseq;
    } else if ( Brev ) {
      B = Brev.get();
    } else {
      Brev.reset(new char[Blen + 2]);
      //      strcpy ( Brev.get() + 1, Bseq + 1 );
      memcpy ( Brev.get() + 1, Bseq + 1, Blen);
      Brev[0] = Brev[Blen + 1] = '\0';
      Reverse_Complement (Brev.get(), 1, Blen);
      B = Brev.get();
    }

    for ( Mp = Cp->matches.begin( ); Mp < Cp->matches.end( ); ++Mp) {
      //-- always_assert for each match in cluster, it is indeed a match
      long int x = Mp->sA;
      long int y = Mp->sB;
      for (long int i = 0; i < Mp->len; i ++ )
        always_assert ( A[x ++] == B[y ++] );

      //-- always_assert for each match in cluster, it is contained in an alignment
      for ( Ap = Alignments.begin( ); Ap < Alignments.end( ); Ap ++ ) {
        if ( Ap->sA <= Mp->sA  &&  Ap->sB <= Mp->sB  &&
             Ap->eA >= Mp->sA + Mp->len - 1  &&
             Ap->eB >= Mp->sB + Mp->len - 1 )
          break;
      }
      always_assert ( Ap < Alignments.end( ) );
    }
  }

  //-- always_assert alignments are optimal (quick check if first and last chars equal)
  for ( Ap = Alignments.begin( ); Ap < Alignments.end( ); ++Ap) {
    if ( Ap->dirB == REVERSE_CHAR ) {
      always_assert (Brev);
      B = Brev.get();
    } else
      B = Bseq;
    always_assert ( Ap->sA <= Ap->eA );
    always_assert ( Ap->sB <= Ap->eB );

    always_assert ( Ap->sA >= 1 && Ap->sA <= Alen );
    always_assert ( Ap->eA >= 1 && Ap->eA <= Alen );
    always_assert ( Ap->sB >= 1 && Ap->sB <= Blen );
    always_assert ( Ap->eB >= 1 && Ap->eB <= Blen );

    char Xc = toupper(isalpha(A[Ap->sA]) ? A[Ap->sA] : sw_align::STOP_CHAR);
    char Yc = toupper(isalpha(B[Ap->sB]) ? B[Ap->sB] : sw_align::STOP_CHAR);
    always_assert ( 0 <= sw_align::MATCH_SCORE [0] [Xc - 'A'] [Yc - 'A'] );

    Xc = toupper(isalpha(A[Ap->eA]) ? A[Ap->eA] : sw_align::STOP_CHAR);
    Yc = toupper(isalpha(B[Ap->eB]) ? B[Ap->eB] : sw_align::STOP_CHAR);
    always_assert ( 0 <= sw_align::MATCH_SCORE [0] [Xc - 'A'] [Yc - 'A'] );
  }
}

void printDeltaAlignments(const std::vector<Alignment>& Alignments,
                          const std::string& AId, const long Alen,
                          const std::string& BId, const long Blen,
                          std::ostream& DeltaFile, const long minLen)
//  Simply output the delta information stored in Alignments to the
//  given delta file. Free the memory used by Alignments once the
//  data is successfully output to the file.
{
  bool header = false;
  for(const auto& A : Alignments) {
    if(std::abs(A.eA - A.sA) + 1 < minLen && std::abs(A.eB - A.sB) + 1 < minLen)
      continue;
    if(!header) {
      DeltaFile << '>' << AId << ' ' << BId << ' ' << Alen << ' ' << Blen << '\n';
      header = true;
    }
    const bool fwd = A.dirB == FORWARD_CHAR;
    DeltaFile << A.sA << ' ' << A.eA << ' '
              << (fwd ? A.sB : revC(A.sB, Blen)) << ' '
              << (fwd ? A.eB : revC(A.eB, Blen)) << ' '
              << A.Errors << ' ' << A.SimErrors << ' ' << A.NonAlphas
              << '\n';

    for(const auto& D : A.delta)
      DeltaFile << D << '\n';
    DeltaFile << "0\n";
  }
}

std::string createCIGAR(const std::vector<long int>& ds, long int start, long int end, long int len,
                        bool hard_clip) {
  std::string res;
  long int    off   = 0;
  long int    range = 0;
  if(start > 1) {
    res += std::to_string(start - 1) + (hard_clip ? 'H' : 'S');
    off += start - 1;
  }
  for(const auto& id : ds) {
    if(std::abs(id) == 1) {
      // Add together multiple insertions and deletions
      if(range == 0 || (id < 0 && range < 0) || (id > 0 && range > 0)) {
        range += id;
        continue;
      }
    }
    if(range) {
      res += std::to_string(std::abs(range)) + (range > 0 ? 'D' : 'I');
      if(range < 0)
        off += std::abs(range);
      range = 0;
    }
    res += std::to_string(std::abs(id) - 1) + 'M';
    off += std::abs(id) - 1;
    range = (id > 0 ? 1 : -1);
    assert(off <= end);
  }
  if(range) {
    res += std::to_string(std::abs(range)) + (range > 0 ? 'D' : 'I');
    if(range < 0)
      off += std::abs(range);
  }
  if(off < end)
    res += std::to_string(end - off) + 'M';
  if(end < len)
    res += std::to_string(len - end) + (hard_clip ? 'H' : 'S');
  return res;
}

// Create MD string for SAM format.
std::string createMD(const Alignment& al, const char* ref,
                     const char* qry, size_t qry_len) {
  auto        it          = error_iterator_type(al, ref, qry, qry_len);
  const auto  it_end      = error_iterator_type(al, ref);
  bool        in_deletion = false;
  long        pos         = 0;
  long        prev_dst    = 0;
  std::string res;

  for( ; it != it_end; ++it) {
    const long diff = it->dst - prev_dst;
    switch(it->type) {
    case NONE: break; // Error! Should not happen! Ignore for now.
    case MISMATCH:
      res          += std::to_string(diff - 1) + *it->ref;
      in_deletion  = false;
      pos += diff;
      prev_dst = it->dst;
      break;
    case INSERTION:
      prev_dst = 0;
      if(!in_deletion || it->dst > 1) {
        res += std::to_string(diff - 1) + '^' + *it->ref;
        in_deletion = true;
      } else {
        res += *it->ref;
      }
      pos += diff;
      break;
    case DELETION:
      // Kind of ignore that event. Nothing is reported in the MD
      // string, it looks like a match.
      prev_dst     = -diff + 1;
      in_deletion  = false;
      break;
    }
  }
  // if(end < pos) error!
  res += std::to_string(al.eA - pos);

  return res;
}

} // namespace postnuc
} // namespace mummer
