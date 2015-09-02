#include <iostream>
#include <fstream>
#include <climits>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <mummer/nucmer.hpp>
#include <src/umd/simplenucmer_cmdline.hpp>

struct getrealpath {
  const char *path, *res;
  getrealpath(const char* p) : path(p), res(realpath(p, nullptr)) { }
  ~getrealpath() { free((void*)res); }
  operator const char*() const { return res ? res : path; }
};

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
  while(std::getline(is, line)) {
    for(const char base : line)
      res += std::tolower(base);
    //    res += line;
  }
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
  getrealpath real_ref(args.ref_arg), real_qry(args.qry_arg);
  os << real_ref << ' ' << real_qry << '\n'
     << "NUCMER\n";

  std::mutex os_mutex;
  mummer::nucmer::FileAligner aligner(args.ref_arg, opts);
  auto print_function = [&](std::vector<mummer::postnuc::Alignment>&& als,
                            const mummer::nucmer::FastaRecordPtr& Af, const mummer::nucmer::FastaRecordSeq& Bf) {
    std::lock_guard<std::mutex> lock (os_mutex);
    mummer::postnuc::printDeltaAlignments(als, Af.Id(), Af.len(), Bf.Id(), Bf.len(), os);
    if(!os.good())
      simplenucmer_cmdline::error() << "Error while writing to output delta file '" << args.delta_arg << '\'';
  };
  os << std::flush;

  const unsigned int nb_threads = args.threads_given ? args.threads_arg : std::thread::hardware_concurrency();
  aligner.align_file(args.qry_arg, print_function, nb_threads);
  return 0;
}
