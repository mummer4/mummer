#include <iostream>
#include <fstream>
#include <climits>
#include <cstdlib>
#include <thread>
#include <memory>
#include <mummer/nucmer.hpp>
#include <src/umd/nucmer_cmdline.hpp>
#include <thread_pipe.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

struct getrealpath {
  const char *path, *res;
  getrealpath(const char* p) : path(p), res(realpath(p, nullptr)) { }
  ~getrealpath() { free((void*)res); }
  operator const char*() const { return res ? res : path; }
};

typedef std::vector<const char*>::const_iterator         path_iterator;
typedef jellyfish::stream_manager<path_iterator>         stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;

void query_thread(mummer::nucmer::FileAligner* aligner, sequence_parser* parser,
                  thread_pipe::ostream_buffered* printer, const nucmer_cmdline* args) {
  auto output_it = printer->begin();
  const bool sam = args->sam_short_given || args->sam_long_given;

  auto print_function = [&](std::vector<mummer::postnuc::Alignment>&& als,
                            const mummer::nucmer::FastaRecordPtr& Af, const mummer::nucmer::FastaRecordSeq& Bf) {
    assert(Af.Id()[strlen(Af.Id()) - 1] != ' ');
    assert(Bf.Id().back() != ' ');
    if(!sam)
      mummer::postnuc::printDeltaAlignments(als, Af.Id(), Af.len(), Bf.Id(), Bf.len(), *output_it, args->minalign_arg);
    else
      mummer::postnuc::printSAMAlignments(als, Af, Bf, *output_it, args->sam_long_given, args->minalign_arg);
    if(output_it->tellp() > 1024)
      ++output_it;
  };
  aligner->thread_align_file(*parser, print_function);
  output_it.done();
}

void query_long(mummer::nucmer::FileAligner* aligner, sequence_parser* parser,
                thread_pipe::ostream_buffered* printer, const nucmer_cmdline* args) {
  auto output_it = printer->begin();
  auto print_function = [&](std::vector<mummer::postnuc::Alignment>&& als,
                            const mummer::nucmer::FastaRecordPtr& Af, const mummer::nucmer::FastaRecordSeq& Bf) {
    mummer::postnuc::printDeltaAlignments(als, Af.Id(), Af.len(), Bf.Id(), Bf.len(), *output_it, args->minalign_arg);
    if(output_it->tellp() > 1024)
      ++output_it;
  };

  while(true) {
    sequence_parser::job j(*parser);
    if(j.is_empty()) break;
    for(size_t i = 0; i < j->nb_filled; ++i) {
      mummer::nucmer::FastaRecordSeq Query(j->data[i].seq.c_str(), j->data[i].seq.length(), j->data[i].header.c_str());
      aligner->align_long_sequences(Query, print_function);
    }
  }
  output_it.done();
}

int main(int argc, char *argv[]) {
  std::ios::sync_with_stdio(false);
  std::string cmdline(argv[0]); // Save command line
  for(int i = 1; i < argc; ++i) {
    cmdline += " ";
    cmdline += argv[i];
  }

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

  const std::string output_file =
    args.delta_given ? args.delta_arg
    : (args.sam_short_given ? args.sam_short_arg
       : (args.sam_long_given ? args.sam_long_arg
          : args.prefix_arg + ".delta"));
  std::ofstream os;
  if(!args.qry_arg.empty()) {
    if(args.qry_arg.size() != 1 && !(args.sam_short_given || args.sam_long_given))
      nucmer_cmdline::error() << "Multiple query file is only supported with the SAM output format";
    os.open(output_file);
    if(!os.good())
      nucmer_cmdline::error() << "Failed to open output file '" << output_file << '\'';

    getrealpath real_ref(args.ref_arg), real_qry(args.qry_arg[0]);
    if(args.sam_short_given || args.sam_long_given) {
      os << "@HD VN1.0 SO:unsorted\n"
         << "@PG ID:nucmer PN:nucmer VN:4.0 CL:\"" << cmdline << "\"\n";
    } else {
      os << real_ref << ' ' << real_qry << '\n'
         << "NUCMER\n";
    }
  }
  thread_pipe::ostream_buffered output(os);

  std::unique_ptr<mummer::nucmer::FileAligner> aligner;
  std::ifstream reference;

  if(args.load_given) {
    mummer::nucmer::sequence_info reference_info(args.ref_arg);
    mummer::mummer::sparseSA SA(reference_info.sequence, args.load_arg);
    aligner.reset(new mummer::nucmer::FileAligner(std::move(reference_info), std::move(SA), opts));
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

    stream_manager     streams(args.qry_arg.cbegin(), args.qry_arg.cend());
    const unsigned int nb_threads = args.threads_given ? args.threads_arg : std::thread::hardware_concurrency();
#ifdef _OPENMP
    if(args.threads_given) omp_set_num_threads(nb_threads);
#endif // _OPENMP

    if(!args.genome_flag) {
      sequence_parser    parser(4 * nb_threads, 10, args.max_chunk_arg, 1, streams);

#ifdef _OPENMP
#pragma omp parallel
      {
        query_thread(aligner.get(), &parser, &output, &args);
      }
#else // _OPENMP
      std::vector<std::thread> threads;
      for(unsigned int i = 0; i < nb_threads; ++i)
        threads.push_back(std::thread(query_thread, aligner.get(), &parser, &output, &args));

      for(auto& th : threads)
        th.join();
#endif // _OPENMP
    } else {
      // Genome flag on
      sequence_parser    parser(4, 1, 1, streams);
      query_long(aligner.get(), &parser, &output, &args);
    }
  } while(!args.load_given && reference.peek() != EOF);
  output.close();
  os.close();

  return 0;
}
