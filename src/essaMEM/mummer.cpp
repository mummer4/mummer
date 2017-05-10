#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <mummer/sparseSA.hpp>
#include <mummer/fasta.hpp>
#include <thread_pipe.hpp>

#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <cctype> // std::tolower(), uppercase/lowercase conversion

// NOTE use of special characters ~, `, and $ !!!!!!!!

// To read input in parallel
typedef jellyfish::stream_manager<const char**>          stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;


void usage(std::string prog);

enum mum_t { MUM, MAM, MEM };

int   min_len          = 20;
int   sparseMult       = 1;
mum_t type             = MAM;
bool  rev_comp         = false;
bool  _4column         = false;
bool  nucleotides_only = false;
bool  forward          = true;
bool  setRevComp       = false;
bool  setBoth          = false;
bool  automatic        = true;
bool  automaticSkip    = true;
bool  automaticKmer    = true;
bool  suflink          = true;
bool  child            = false;
int   kmer             = 0;
bool  print_length     = false;
bool  printSubstring   = false;
bool  printRevCompForw = false;
int   K                = 1;
int   num_threads      = 1;
int   query_threads    = 0;
size_t max_chunk       = 50000;

// Information on matches and write in multi-thread safe
struct match_info {
  std::string                          meta;
  size_t                               len;
  std::vector<mummer::mummer::match_t> fwd_matches, bwd_matches;
  void clear() {
    meta.clear();
    fwd_matches.clear();
    bwd_matches.clear();
  }
};

void print_match_info(std::ostream& os, const match_info& m, const mummer::mummer::sparseSAMatch* sa) {
    os << "> " << m.meta;
    if(print_length) os << "\tLen = " << m.len;
    os << '\n';
    for(const auto& i : m.fwd_matches)
      sa->print_match(os, i);
    if(rev_comp) {
      os << "> " << m.meta << " Reverse";
      if(print_length) os << "\tLen = " << m.len;
      os << '\n';
      for(const auto& i : m.bwd_matches)
        sa->print_match(os, i);
    }
}

void query_thread(const mummer::mummer::sparseSAMatch* sa, sequence_parser* parser,
                  thread_pipe::ostream_buffered* printer) {
  auto       output_it = printer->begin();
  match_info match;

  while(true) {
    // Get a job (a batch of sequences)
    sequence_parser::job j(*parser);
    if(j.is_empty()) break;

    // Process each sequence in job
    for(size_t i = 0; i < j->nb_filled; ++i, ++output_it) {

      // Get meta (header)
      match.clear();
      std::string& header = j->data[i].header;
      size_t start_meta = header.find_first_not_of(" ");
      if(start_meta == std::string::npos) continue;
      size_t end_meta = header.find_first_of(" ", start_meta);
      match.meta = header.substr(start_meta, std::min(end_meta, header.size()) - start_meta);

      // Clean up sequence if nucleotides only
      std::string& P = j->data[i].seq;
      match.len = P.size();
      char* const end = (char*)P.data() + P.size();
      if(nucleotides_only) {
        for(char* seq = (char*)P.data(); seq != end; ++seq) {
          switch(*seq) {
          case 'a': case 't': case 'g': case 'c': break;
          case 'A': *seq = 'a'; break;
          case 'C': *seq = 'c'; break;
          case 'G': *seq = 'g'; break;
          case 'T': *seq = 't'; break;
          default:
            *seq = '~';
          }
        }
      } else { // !nucleotides_only
        for(char* seq = (char*)P.data(); seq != end; ++seq)
          *seq = std::tolower(*seq);
      }

      // Get matches
      if(forward) {
        match.fwd_matches.clear();
        switch(type) {
        case MAM: sa->MAM(P, min_len, false, match.fwd_matches); break;
        case MUM: sa->MUM(P, min_len, false, match.fwd_matches); break;
        case MEM: sa->MEM(P, min_len, false, match.fwd_matches); break;
        }
      }
      if(rev_comp) {
        match.bwd_matches.clear();
        reverse_complement(P, nucleotides_only);
        switch(type) {
        case MAM: sa->MAM(P, min_len, printRevCompForw, match.bwd_matches); break;
        case MUM: sa->MUM(P, min_len, printRevCompForw, match.bwd_matches); break;
        case MEM: sa->MEM(P, min_len, printRevCompForw, match.bwd_matches); break;
        }
      }
      print_match_info(*output_it, match, sa);
    }
  }
  output_it.done();
}

int main(int argc, char* argv[]) {
  std::ios::sync_with_stdio(false);

  // Collect parameters from the command line.
  std::string save;
  std::string load;

  while (1) {
    static struct option long_options[] = {
      {"l", 1, 0, 0}, // 0
      {"mumreference", 0, 0, 0}, // 1
      {"b", 0, 0, 0}, // 2
      {"maxmatch", 0, 0, 0}, // 3
      {"mum", 0, 0, 0}, // 4
      {"mumcand", 0, 0, 0},  // 5
      {"F", 0, 0, 0}, // 6
      {"k", 1, 0, 0}, // 7
      {"threads", 1, 0, 0}, // 8
      {"n", 0, 0, 0}, // 9
      {"qthreads", 1, 0, 0}, // 10
      {"suflink", 1, 0, 0}, // 11
      {"child", 1, 0, 0}, // 12
      {"skip", 1, 0, 0}, // 13
      {"L", 0, 0, 0}, // 14
      {"r", 0, 0, 0}, // 15
      {"s", 0, 0, 0}, // 16
      {"c", 0, 0, 0}, // 17
      {"kmer", 1, 0, 0}, // 18
      {"save", 1, 0, 0}, // 19
      {"load", 1, 0, 0}, // 20
      {"max-chunk", 1, 0, 0}, // 21
      {"version", 0, 0, 0}, // 22
      {0, 0, 0, 0}
    };
    int longindex = -1;
    int c = getopt_long_only(argc, argv, "", long_options, &longindex);
    if(c == -1) break; // Done parsing flags.
    else if(c == '?') { // If the user entered junk, let him know.
      std::cerr << "Invalid parameters." << std::endl;
      usage(argv[0]);
    }
    else {
      // Branch on long options.
      switch(longindex) {
      case 0: min_len = atol(optarg); break;
      case 1: type = MAM; break;
      case 2: setBoth = true;	break;
      case 3: type = MEM; break;
      case 4: type = MUM; break;
      case 5: type = MAM; break;
      case 6: _4column = true; break;
      case 7: K = atoi(optarg); break;
      case 8: num_threads = atoi(optarg); break;
      case 9: nucleotides_only = true; break;
      case 10: query_threads = atoi(optarg) ; break;
      case 11: suflink = atoi(optarg) > 0;	automatic = false; break;
      case 12: child = atoi(optarg) > 0;	automatic = false; break;
      case 13: sparseMult = atoi(optarg); automaticSkip = false; break;
      case 14: print_length = true; break;
      case 15: setRevComp = true; break;
      case 16: printSubstring = true; break;
      case 17: printRevCompForw = true; break;
      case 18: kmer = atoi(optarg); automaticKmer = false; break;
      case 19: save = optarg; break;
      case 20: load = optarg; break;
      case 21: max_chunk = atoi(optarg); break;
      case 22:
#ifdef VERSION
        std::cout << VERSION << '\n';
#else
        std::cout << "<unknown version>\n";
#endif
        exit(0);
      default: break;
      }
    }
  }
  if (argc - optind < 2) usage(argv[0]);

  if(K != 1 && type != MEM) { std::cerr << "-k option valid only for -maxmatch" << std::endl; exit(1); }
  if(num_threads <= 0) { std::cerr << "invalid number of threads specified" << std::endl; exit(1); }
  if(query_threads <= 0) { query_threads = std::thread::hardware_concurrency(); }

  std::string ref_fasta = argv[optind];
  int argNumber = optind+1;
  // numQueryFiles = 0;
  // while(argNumber < argc){
  //     query_fasta[numQueryFiles] = argv[argNumber];
  //     numQueryFiles++;
  //     argNumber++;
  // }

  std::string ref;

  std::vector<std::string> refdescr;
  std::vector<long> startpos;

  load_fasta(ref_fasta, ref, refdescr, startpos);

  // Automatically use 4 column format if there are multiple reference sequences.
  if(startpos.size() > 1) _4column = true;
  if(automatic){
      suflink = K < 4;
      child = K >= 4;
  }
  if(automaticSkip){
      if(suflink && !child) sparseMult = 1;
      else{
          if(K >= 4) sparseMult = (int) std::max((min_len-10)/K,1);
          else sparseMult = (int) std::max((min_len-12)/K,1);
      }
  }
  else{
      if(sparseMult*K > min_len){
        while(sparseMult*K > min_len)
            sparseMult--;
        std::cerr << "skip parameter was decreased to " << sparseMult << " because skip*K > minimum length" << std::endl;
      }
      if(sparseMult*K > min_len-10){
          std::cerr << "note that the skip parameter is very high, a value of " << ((int) (min_len-10)/K);
          std::cerr << " or " << ((int) (min_len-12)/K) << " would be more appropriate" << std::endl;
      }
  }
  if(automaticKmer){
      kmer = std::max(0,std::min(10,min_len - sparseMult*K + 1));
  }
  else{
      if(kmer > 12)
          std::cerr << "warning: very large value for kmer-size: index will be very large" << std::endl;
      if(kmer > min_len - sparseMult*K + 1){
          kmer = std::max(0,std::min(10,min_len - sparseMult*K + 1));
        std::cerr << "kmer size was reduced to " << kmer << " because the user set value is too large and cannotbe used in the algorithm" << std::endl;
      }
  }

  if(setBoth && setRevComp){
      std::cerr << "ERROR -r and -b options are mutually exclusive" << std::endl;
      exit(1);
  }
  if(setBoth || setRevComp)
      rev_comp = true;
  if(setRevComp)
      forward = false;
  std::unique_ptr<mummer::mummer::sparseSAMatch> sa(new mummer::mummer::sparseSAMatch(ref, refdescr, startpos, _4column, K, suflink, child, kmer>0, sparseMult, kmer, printSubstring, nucleotides_only));
  if(!load.empty()){
    if(sa->load(load)){
      std::cerr << "index loaded succesfully\n"
                << "WARNING: program does not check the soundness of the reference file for given loaded index. Use the same reference file as used for constructing the index\n"
                << "WARNING: some options are now taken from loaded index instead of current user-set values\n"
                << "these include: sparseness (-k), suffix links (-suflink), child array (-child) and kmer table size (-kmer)." << std::endl;
          //update sparseMult if necessary
          if(automaticSkip){
            if(sa->hasSufLink && !sa->hasChild) sparseMult = 1;
            else{
                if(sa->K >= 4) sparseMult = (int) std::max((min_len-10)/sa->K,1L);
                else sparseMult = (int) std::max((min_len-12)/sa->K,1L);
            }
         }
         else{
            if(sparseMult*sa->K > min_len){
                while(sparseMult*sa->K > min_len)
                    sparseMult--;
                std::cerr << "skip parameter was decreased to " << sparseMult << " because skip*K > minimum length" << std::endl;
            }
            if(sparseMult*sa->K > min_len-10){
              std::cerr << "note that the skip parameter is very high, a value of " << ((int) (min_len-10)/K) << '\n'
                        << " or " << ((int) (min_len-12)/sa->K) << " would be more appropriate" << std::endl;
            }
         }
          sa->sparseMult = sparseMult;
         //update other fields
         suflink = sa->hasSufLink;
         child = sa->hasChild;
         kmer = sa->kMerSize;
         K = sa->K;
      }
      else{
          std::cerr << "unable to load index " << load << '\n'
                    << "construct new index..." << std::endl;
          sa->construct();
      }
  }
  else{
      sa->construct();
  }
  if(!save.empty()){
      sa->save(save);
  }

  // Open input files
  stream_manager  streams((const char**)(argv + argNumber), (const char**)(argv + argc));
  sequence_parser               parser(4 * query_threads, 10, max_chunk, 1, streams);
  thread_pipe::ostream_buffered output(std::cout);

  // Launch query threads
  std::vector<std::thread> threads;
  for(int i = 0; i < query_threads; ++i)
    threads.push_back(std::thread(query_thread, sa.get(), &parser, &output));

  // Wait for all threads to terminate.
  for(auto& th : threads)
    th.join();
  output.close();
}


void usage(std::string prog) {
  std::cerr << "Usage: " << prog << " [options] <reference-file> <query file1> . . . [query file32]" << '\n'
            << "Implemented MUMmer v3 options:" << '\n'
            << "-mum           compute maximal matches that are unique in both sequences" << '\n'
            << "-mumreference  compute maximal matches that are unique in the reference-" << '\n'
            << "               sequence but not necessarily in the query-sequence (default)" << '\n'
            << "-mumcand       same as -mumreference" << '\n'
            << "-maxmatch      compute all maximal matches regardless of their uniqueness" << '\n'
            << "-l             set the minimum length of a match" << '\n'
            << "               if not set, the default value is 20" << '\n'
            << "-b             compute forward and reverse complement matches" << '\n'
            << "-F             force 4 column output format regardless of the number of" << '\n'
            << "               reference sequence inputs"  << '\n'
            << "-n             match only the characters a, c, g, or t" << '\n'
            << "-L             print length of query sequence in header of matches" << '\n'
            << "-r             compute only reverse complement matches" << '\n'
            << "-s             print first 53 characters of the matching substring" << '\n'
            << "-c             Report the query position of a reverse complement match relative to the forward strand of the query sequence" << '\n'
            << '\n'
            << "Additional options:" << '\n'
            << "-k             sampled suffix positions (one by default)" << '\n'
            << "-threads       number of threads to use for -maxmatch, only valid k > 1 " << '\n'
            << "-qthreads      number of threads to use for queries " << '\n'
            << "-suflink       use suffix links (1=yes or 0=no) in the index and during search [auto]" << '\n'
            << "-child         use child table (1=yes or 0=no) in the index and during search [auto]" << '\n'
            << "-skip          sparsify the MEM-finding algorithm even more, performing jumps of skip*k [auto (l-10)/k]" << '\n'
            << "               this is a performance parameter that trade-offs SA traversal with checking of right-maximal MEMs" << '\n'
            << "-kmer          use kmer table containing sa-intervals (speeds up searching first k characters) in the index and during search [int value, auto]" << '\n'
            << "-save (string) save index to file to use again later (string)" << '\n'
            << "-load (string) load index from file" << '\n'
            << '\n'
            << "Example usage:" << '\n'
            << '\n'
            << "./mummer -maxmatch -l 20 -b -n -k 3 -threads 3 ref.fa query.fa" << '\n'
            << "Find all maximal matches on forward and reverse strands" << '\n'
            << "of length 20 or greater, matching only a, c, t, or g." << '\n'
            << "Index every 3rd position in the ref.fa and use 3 threads to find MEMs." << '\n'
            << "Fastest method for one long query sequence." << '\n'
            << '\n'
            << "./mummer -maxmatch -l 20 -b -n -k 3 -qthreads 3 ref.fa query.fa" << '\n'
            << "Same as above, but now use a single thread for every query sequence in" << '\n'
            << "query.fa. Fastest for many small query sequences." << std::endl;
  exit(1);
}
