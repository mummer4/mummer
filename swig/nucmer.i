%{
#include <sstream>
%}

%template(LongVector) ::std::vector<long>;

namespace mummer {

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
};
::std::vector<mummer::postnuc::Alignment> align_sequences(const char* reference, const char* query,
                                                          mummer::nucmer::Options opts = mummer::nucmer::Options());
} // namespace nucmer
} // namespace mummer

