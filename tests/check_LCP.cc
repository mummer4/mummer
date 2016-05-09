#include <stdexcept>
#include <fstream>
#include <iostream>

#undef NDEBUG
#include <mummer/sparseSA.hpp>
#include <mummer/fasta.hpp>
#include <mummer/sparseSA_imp.hpp>

#include "check_LCP_cmdline.hpp"

typedef mummer::mummer::sparseSA_aux aux_type;
typedef mummer::mummer::vector_32_48 SA_type;
typedef mummer::mummer::vec_uchar lcp_type;

int main(int argc, char *argv[]) {
  check_LCP_cmdline args(argc, argv);

  aux_type                 aux_info(args.prefix_arg + ".aux");
  SA_type                  SA(args.prefix_arg + ".sa");
  SA_type                  ISA(args.prefix_arg + ".isa");
  lcp_type                 LCP(SA);
  std::string              sequence;
  std::vector<std::string> descr;
  std::vector<long>        startpos;
  load_fasta(args.sequence_arg, sequence, descr, startpos);
  mummer::sparseSA_imp::computeLCP(LCP, sequence, SA, ISA, aux_info.N, aux_info.K);
  LCP.init();
  return 0;
}
