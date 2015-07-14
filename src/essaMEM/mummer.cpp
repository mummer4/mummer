#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "sparseSA.hpp"
#include "fasta.hpp"

#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <cctype> // std::tolower(), uppercase/lowercase conversion

// NOTE use of special characters ~, `, and $ !!!!!!!!

// using namespace std;

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
int   query_threads    = 1;

mummer::mummer::sparseSA *sa;
std::string query_fasta[32];
int MAX_QUERY_FILES = 32;
int numQueryFiles = 0;

struct query_arg {
  int skip0;
  int skip;
  int queryFile;
};

void *query_thread(void *arg_) {
  query_arg *arg = (query_arg *)arg_;
  std::ifstream data(query_fasta[arg->queryFile].c_str());

  std::vector<mummer::mummer::match_t> matches;

  const bool print = arg->skip == 1;


  if(!data.is_open()) { std::cerr << "unable to open " << query_fasta[arg->queryFile] << std::endl; exit(1); }

  int c = data.peek();
  if(c != '>') {
    std::cerr << "error, first character must be a '>', got '" << (char)c << "'" << std::endl;
    exit(1);
  }

  std::string P, meta, line;
  for(long seq_cnt = 0 ; c != EOF; seq_cnt++) {
    P.clear(); meta.clear();

    // Load metadata
    getline(data, line);
    size_t start = line.find_first_not_of(" ", 1);
    if(start != std::string::npos) {
      size_t end = line.find_first_of(" ", start);
      meta = line.substr(start, std::min(end, line.size()) - start);
    }

    // Load sequence
    for(c = data.peek(); c != EOF && c != '>'; c = data.peek()) {
      getline(data, line); // Load one line at a time.
      long start = 0, end = line.length() - 1;
      trim(line, start,end);
      for(long i = start; i <= end; i++) {
        char c = std::tolower(line[i]);
        if(nucleotides_only) {
          switch(c) {
          case 'a': case 't': case 'g': case 'c': break;
          default:
            c = '~';
          }
        }
        P += c;
      }
    }

    if(meta.empty()) continue;
    if(seq_cnt % arg->skip == arg->skip0) {
      // Process P.
      //   std::cerr << "# P.length()=" << P.length() << std::endl;
      if(forward){
        if(print){
          if(print_length) std::cout << "> " << meta << "\tLen = " << P.length() << '\n';
          else std::cout << "> " << meta << '\n';
        }
        switch(type) {
        case MAM: sa->MAM(P, min_len, true, std::cout); break;
        case MUM: sa->MUM(P, min_len, true, std::cout); break;
        case MEM: sa->MEM(P, min_len, true, std::cout); break;
        }
        if(!print) sa->print_match(std::cout, meta, false);
      }
      if(rev_comp) {
        reverse_complement(P, nucleotides_only);
        if(print){
          if(print_length) std::cout << "> " << meta << " Reverse\tLen = " << P.length() << '\n';
          else std::cout << "> " << meta << " Reverse\n";
        }
        switch(type) {
        case MAM: sa->MAM(P, min_len, false, std::cout); break;
        case MUM: sa->MUM(P, min_len, false, std::cout); break;
        case MEM: sa->MEM(P, min_len, false, std::cout); break;
        }
        if(!print) sa->print_match(std::cout, meta, true);
      }
    }
  }

  //  std::cerr << "number of M(E/A/U)Ms: " << memCounter << std::endl;
  pthread_exit(NULL);
  return 0;
}

// Added by Simon Gog for testing
// void write_lock(int i){
//   std::ofstream lockfile("lock.txt", std::ios_base::trunc);
// 	lockfile<<i<<std::endl;
// 	lockfile.close();
// }

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
      default: break;
      }
    }
  }
  if (argc - optind < 2 || argc - optind >  MAX_QUERY_FILES + 1) usage(argv[0]);

  if(K != 1 && type != MEM) { std::cerr << "-k option valid only for -maxmatch" << std::endl; exit(1); }
  if(num_threads <= 0) { std::cerr << "invalid number of threads specified" << std::endl; exit(1); }

  std::string ref_fasta = argv[optind];
  int argNumber = optind+1;
  numQueryFiles = 0;
  while(argNumber < argc){
      query_fasta[numQueryFiles] = argv[argNumber];
      numQueryFiles++;
      argNumber++;
  }

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
 
  sa = new mummer::mummer::sparseSA(ref, refdescr, startpos, _4column, K, suflink, child, kmer>0, sparseMult, kmer, printSubstring, printRevCompForw, nucleotides_only);
  if(!load.empty()){
      std::cerr << "attempting to load index " << load << std::endl;
      if(sa->load(load)){
          std::cerr << "index loaded succesfully" << std::endl;
          std::cerr << "WARNING: program does not check the soundness of the reference file for given loaded index. Use the same reference file as used for constructing the index" << std::endl;
          std::cerr << "WARNING: some options are now taken from loaded index instead of current user-set values" << std::endl;
          std::cerr << "these include: sparseness (-k), suffix links (-suflink), child array (-child) and kmer table size (-kmer)." << std::endl;
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
                std::cerr << "note that the skip parameter is very high, a value of " << ((int) (min_len-10)/K);
                std::cerr << " or " << ((int) (min_len-12)/sa->K) << " would be more appropriate" << std::endl;
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
          std::cerr << "unable to load index " << load << std::endl;
          std::cerr << "construct new index..." << std::endl;
          sa->construct();
      }
  }
  else{
      sa->construct();
  }
  if(!save.empty()){
      sa->save(save);
      std::cerr << "saved index to " << save << std::endl;
  }
  std::cerr << "INDEX SIZE IN BYTES: " << sa->index_size_in_bytes() << std::endl;

  clock_t start = clock();
  rusage m_ruse1, m_ruse2;
  getrusage(RUSAGE_SELF, &m_ruse1);
  pthread_attr_t attr;  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  for(int idx = 0; idx < numQueryFiles; idx++){
    std::vector<query_arg> args(query_threads);
    std::vector<pthread_t> thread_ids(query_threads); 

    // Initialize additional thread data.
    for(int i = 0; i < query_threads; i++) {
        args[i].skip = query_threads;
        args[i].skip0 = i;
        args[i].queryFile = idx;
    }

    // Create joinable threads to find MEMs.
    for(int i = 0; i < query_threads; i++)
        pthread_create(&thread_ids[i], &attr, query_thread, (void *)&args[i]);

    // Wait for all threads to terminate.
    for(int i = 0; i < query_threads; i++)
        pthread_join(thread_ids[i], NULL);   
  }
  clock_t end = clock();
  getrusage(RUSAGE_SELF, &m_ruse2);
  double wall_time = (double)( end - start ) /CLOCKS_PER_SEC;
  std::cerr << "mapping: done" << std::endl;
  std::cerr << "time for mapping (wall time): " << wall_time << std::endl;
  timeval t1, t2;
  t1 = m_ruse1.ru_utime;
  t2 = m_ruse2.ru_utime;
  double cpu_time = ((double)(t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec )))/1000.0;
  std::cerr << "time for mapping (cpu time): " << cpu_time << std::endl;
  t1 = m_ruse1.ru_stime;
  t2 = m_ruse2.ru_stime;
  double sys_time = ((double)(t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec )))/1000.0;
  std::cerr << "time for mapping (sys time): " << sys_time << std::endl;

  delete sa;
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
