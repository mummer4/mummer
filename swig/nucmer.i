%{
#include <sstream>
%}

%template(LongVector) ::std::vector<long>;

namespace mummer {

void set_num_threads(int nb);
int  get_num_threads();

namespace postnuc {
struct Alignment   {
  //-- An alignment object between two sequences A and B
  signed char         dirB;     // the query sequence direction
  long int            sA, sB, eA, eB; // the start in A, B and the end in A, B
  ::std::vector<long> delta;    // the delta values, with NO zero at the end
  long int            deltaApos; // sum of abs(deltas) - #of negative deltas
  long int            Errors, SimErrors, NonAlphas; // errors, similarity errors, nonalphas
  double identity() const;
  double similarity() const;
  double stopity() const;
  %extend {
    %feature("autodoc", "Return string representation of alignment");
    std::string __str__() {
      std::ostringstream os;
      os << *$self;
      return os.str();
    }
  }
};
%template(AlignmentVector) ::std::vector<mummer::postnuc::Alignment>;
} // namespace postnuc

namespace nucmer {
struct Options {
  mummer::nucmer::Options& mum();
  mummer::nucmer::Options& mumcand();
  mummer::nucmer::Options& mumreference();
  mummer::nucmer::Options& maxmatch();
  mummer::nucmer::Options& breaklen(long l);
  mummer::nucmer::Options& banded();
  mummer::nucmer::Options& nobanded();
  mummer::nucmer::Options& mincluster(long m);
  mummer::nucmer::Options& diagdiff(long d);
  mummer::nucmer::Options& diagfactor(double f);
  mummer::nucmer::Options& extend();
  mummer::nucmer::Options& noextend();
  mummer::nucmer::Options& forward();
  mummer::nucmer::Options& maxgap(long m);
  mummer::nucmer::Options& minmatch(long m);
  mummer::nucmer::Options& optimize();
  mummer::nucmer::Options& nooptimize();
  mummer::nucmer::Options& reverse();
  mummer::nucmer::Options& simplify();
  mummer::nucmer::Options& nosimplify();

  // Options for mummer
  //  match_type match;
  int        min_len;
  //  ori_type   orientation;

  // Options for mgaps
  long   fixed_separation;
  long   max_separation;
  long   min_output_score;
  double separation_factor;
  bool   use_extent;

  // Options for postnuc
  bool do_delta;
  bool do_extend;
  bool to_seqend;
  bool do_shadows;
  int  break_len;
  int  banding;
};
%apply (const char* STRING, size_t LENGTH) { (const char* reference, size_t reference_len) };
%apply (const char* STRING, size_t LENGTH) { (const char* query, size_t query_len) };
::std::vector<mummer::postnuc::Alignment> align_sequences(const char* reference, size_t reference_len,
                                                          const char* query, size_t query_len,
                                                          mummer::nucmer::Options opts = mummer::nucmer::Options());
} // namespace nucmer
} // namespace mummer

