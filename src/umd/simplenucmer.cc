#include <iostream>
#include <fstream>
#include <mummer/nucmer.hpp>
#include <src/umd/simplenucmer_cmdline.hpp>

std::string read_sequence(const char* file, std::string& header) {
  std::string res;
  std::ifstream is(file);
  if(!is.good())
    simplenucmer_cmdline::error() << "Failed to open file '" << file << '\'';
  std::string line;

  std::getline(is, line);
  if(line.size() > 1)
    header = line.substr(1);
  else
    header.clear();
  while(std::getline(is, line))
    res += line;
  return res;
}

int main(int argc, char *argv[]) {
  simplenucmer_cmdline args(argc, argv);
  mummer::nucmer::Options opts;
  opts.breaklen(args.breaklen_arg)
    .mincluster(args.mincluster_arg)
    .diagdiff(args.diagdiff_arg)
    .diagfactor(args.diagfactor_arg)
    .maxgap(args.maxgap_arg)
    .minmatch(args.minmatch_arg);
  if(args.noextend_flag) opts.noextend();
  if(args.nooptimize_flag) opts.nooptimize();
  if(args.nosimplify_flag) opts.nosimplify();
  if(args.forward_flag) opts.forward();
  if(args.reverse_flag) opts.reverse();

  std::ostream os(std::cout.rdbuf());
  std::ofstream delta;
  if(args.delta_given) {
    delta.open(args.delta_arg);
    if(!delta.good())
      simplenucmer_cmdline::error() << "Failed to open output delta file '" << args.delta_arg << '\'';
    os.rdbuf(delta.rdbuf());
  }

  std::string qry_header;
  std::string qry = read_sequence(args.qry_arg, qry_header);

  mummer::nucmer::FileAligner aligner(args.ref_arg, opts);
  aligner.align(qry,
                [&](std::vector<mummer::postnuc::Alignment>&& als,
                   const mummer::nucmer::FastaRecordPtr& Af, const mummer::nucmer::FastaRecordSeq & Bf) {
                  //  const auto alignments = mummer::nucmer::align_sequences(ref.c_str(), qry.c_str());
                  mummer::postnuc::printDeltaAlignments(als, Af.Id(), Af.len(), qry_header, qry.size(),
                                                        os);
                  if(!os.good())
                    simplenucmer_cmdline::error() << "Failed to open output delta file '" << args.delta_arg << '\'';
                });
  return 0;
}
