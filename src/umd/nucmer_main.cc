#include <iostream>
#include <fstream>
#include <climits>
#include <cstdlib>
#include <thread>
#include <memory>
#include <mummer/nucmer.hpp>
#include <src/umd/nucmer_cmdline.hpp>
#include <thread_pipe.hpp>

struct getrealpath {
  const char *path, *res;
  getrealpath(const char* p) : path(p), res(realpath(p, nullptr)) { }
  ~getrealpath() { free((void*)res); }
  operator const char*() const { return res ? res : path; }
};

typedef jellyfish::stream_manager<const char**>          stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;

void query_thread(mummer::nucmer::FileAligner* aligner, sequence_parser* parser,
                  thread_pipe::ostream_buffered* printer, uint32_t minalign) {
  auto output_it = printer->begin();
  auto print_function = [&](std::vector<mummer::postnuc::Alignment>&& als,
                            const mummer::nucmer::FastaRecordPtr& Af, const mummer::nucmer::FastaRecordSeq& Bf) {
    assert(Af.Id()[strlen(Af.Id()) - 1] != ' ');
    assert(Bf.Id().back() != ' ');
    mummer::postnuc::printDeltaAlignments(als, Af.Id(), Af.len(), Bf.Id(), Bf.len(), *output_it, minalign);
    if(output_it->tellp() > 1024)
      ++output_it;
  };
  aligner->thread_align_file(*parser, print_function);
  output_it.done();
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
  if(args.mum_flag) opts.mum();
  if(args.maxmatch_flag) opts.maxmatch();

  const std::string delta_file = args.delta_given ? args.delta_arg : args.prefix_arg + ".delta";
  std::ofstream os(delta_file);
  if(!os.good())
    nucmer_cmdline::error() << "Failed to open output delta file '" << delta_file << '\'';

  getrealpath real_ref(args.ref_arg), real_qry(args.qry_arg);
  os << real_ref << ' ' << real_qry << '\n'
     << "NUCMER\n";
  thread_pipe::ostream_buffered output(os);

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
  do {
    if(!args.load_given)
      aligner.reset(new mummer::nucmer::FileAligner(reference, batch_size,  opts));

    if(args.save_given && !aligner->sa().save(args.save_arg))
      nucmer_cmdline::error() << "Can't save the suffix array to '" << args.save_arg << "'";

    //    os << std::flush;

    const unsigned int nb_threads = args.threads_given ? args.threads_arg : std::thread::hardware_concurrency();
    stream_manager     streams(&args.qry_arg, &args.qry_arg + 1);
    sequence_parser    parser(4 * nb_threads, 10, 1, streams);

    std::vector<std::thread> threads;
    for(unsigned int i = 0; i < nb_threads; ++i)
      threads.push_back(std::thread(query_thread, aligner.get(), &parser, &output, args.minalign_arg));

    for(auto& th : threads)
      th.join();
  } while(!args.load_given && reference.peek() != EOF);
  output.close();
  os.close();

  return 0;
}
