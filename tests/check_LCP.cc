#include <stdexcept>
#include <fstream>
#include <iostream>
#include <utility>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#undef NDEBUG
#include <mummer/sparseSA.hpp>
#include <mummer/sparseSA_imp.hpp>
#include <mummer/nucmer.hpp>

#include "check_LCP_cmdline.hpp"

typedef mummer::mummer::sparseSA_aux   aux_type;
typedef mummer::mummer::vector_32_48   SA_type;
typedef mummer::mummer::vec_uchar      lcp_type;
typedef mummer::nucmer::sequence_info  sequence_type;
typedef mummer::mummer::bounded_string bounded_type;

template<typename T>
void check_LCP_size(const std::string& prefix) {
  T sizeLCP, sizeM;
  { std::ifstream is(prefix + ".lcp");
    is.read((char*)&sizeLCP, sizeof(sizeLCP));
    is.read((char*)&sizeM, sizeof(sizeM));
  }
  struct stat lcp_stat;
  if(stat((prefix + ".lcp").c_str(), &lcp_stat) == -1)
    check_LCP_cmdline::error() << "Failed to stat lcp file";
  const size_t expected = sizeof(sizeLCP) + sizeof(sizeM) + sizeLCP + sizeM * sizeof(mummer::mummer::vec_uchar::item_t);
  if(expected != (size_t)lcp_stat.st_size)
    check_LCP_cmdline::error() << "LCP file is " << lcp_stat.st_size
                               << " but expected " << expected;
}

int main(int argc, char *argv[]) {
  check_LCP_cmdline args(argc, argv);

  if(args.int_flag)
    check_LCP_size<unsigned int>(args.prefix_arg);
  else
    check_LCP_size<size_t>(args.prefix_arg);

  if(args.short_flag)
    return 0;

  aux_type      aux_info(args.prefix_arg + ".aux");
  SA_type       SA(args.prefix_arg + ".sa");
  SA_type       ISA(args.prefix_arg + ".isa");
  lcp_type      LCP(SA);
  lcp_type      LCP_load(SA);
  sequence_type sequence(args.sequence_arg.c_str());
  bounded_type  bounded(sequence.sequence, 1);

  if(aux_info.N <= 0)
    check_LCP_cmdline::error() << "Got invalid N" << aux_info.N;

  if(bounded.size() != (size_t)aux_info.N)
    check_LCP_cmdline::error() << "Sequence len " << sequence.sequence.size() << " != N " << aux_info.N;
  mummer::sparseSA_imp::computeLCP(LCP, bounded, SA, ISA, aux_info.N, aux_info.K);
  LCP.init();
  if(args.output_given) {
    if(!LCP.save(args.output_arg))
      check_LCP_cmdline::error() << "Failed to write LCP";
  }
  if(!LCP_load.load(args.prefix_arg + ".lcp"))
    check_LCP_cmdline::error() << "Failed to load LCP";

  if(LCP_load.vec.size() != LCP.vec.size())
    check_LCP_cmdline::error() << "Vec size differ: " << LCP.vec.size()
                               << " != " << LCP_load.vec.size();
  if(LCP_load.M.size() != LCP.M.size())
    check_LCP_cmdline::error() << "M size differ: " << LCP.M.size()
                               << " != " << LCP_load.M.size();

  if(!std::equal(LCP_load.vec.cbegin(), LCP_load.vec.cend(), LCP.vec.cbegin()))
    check_LCP_cmdline::error() << "Vec differ";

  if(!std::equal(LCP_load.M.cbegin(), LCP_load.M.cend(), LCP.M.cbegin()))
    check_LCP_cmdline::error() << "M differ";

  return 0;
}
