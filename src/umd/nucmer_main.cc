#include <iostream>
#include <fstream>
#include <climits>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <memory>
#include <mummer/nucmer.hpp>
#include <src/umd/nucmer_cmdline.hpp>

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
    nucmer_cmdline::error() << "Failed to open file '" << file << '\'';
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
  nucmer_cmdline args(argc, argv);
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

  const std::string delta_file = args.delta_given ? args.delta_arg : args.prefix_arg + ".delta";
  std::ofstream os(delta_file);
  if(!os.good())
    nucmer_cmdline::error() << "Failed to open output delta file '" << delta_file << '\'';

  getrealpath real_ref(args.ref_arg), real_qry(args.qry_arg);
  os << real_ref << ' ' << real_qry << '\n'
     << "NUCMER\n";

  std::mutex os_mutex;
  std::unique_ptr<mummer::nucmer::FileAligner> aligner;
  std::ifstream reference;

  if(args.load_given) {
    mummer::nucmer::sequence_info reference_info(args.ref_arg);
    mummer::mummer::sparseSA SA(reference_info.sequence, args.load_arg);
    aligner.reset(new mummer::nucmer::FileAligner(std::move(reference_info), std::move(SA)));
  } else {
    reference.open(args.ref_arg);
    if(!reference.good())
      nucmer_cmdline::error() << "Failed to open reference file '" << args.ref_arg << "'";
  }

  const size_t batch_size = args.batch_given ? args.batch_arg : std::numeric_limits<size_t>::max();
  int i = 0;
  do {
    if(!args.load_given)
      aligner.reset(new mummer::nucmer::FileAligner(reference, batch_size,  opts));

    if(args.save_given && !aligner->sa().save(args.save_arg))
      nucmer_cmdline::error() << "Can't save the suffix array to '" << args.save_arg << "'";

    auto print_function = [&](std::vector<mummer::postnuc::Alignment>&& als,
                              const mummer::nucmer::FastaRecordPtr& Af, const mummer::nucmer::FastaRecordSeq& Bf) {
      std::lock_guard<std::mutex> lock (os_mutex);
      assert(Af.Id()[strlen(Af.Id()) - 1] != ' ');
      assert(Bf.Id().back() != ' ');
      mummer::postnuc::printDeltaAlignments(als, Af.Id(), Af.len(), Bf.Id(), Bf.len(), os);
      if(!os.good())
        nucmer_cmdline::error() << "Error while writing to output delta file '" << args.delta_arg << '\'';
    };
    os << std::flush;

    const unsigned int nb_threads = args.threads_given ? args.threads_arg : std::thread::hardware_concurrency();
    aligner->align_file(args.qry_arg, print_function, nb_threads);
    std::cerr << i++ << ' ' << (char)reference.peek() << ' ' << reference.tellg() << std::endl;
  } while(!args.load_given && reference.peek() != EOF);

  return 0;
}
