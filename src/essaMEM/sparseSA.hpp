#ifndef __sparseSA_hpp__
#define __sparseSA_hpp__

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>
#include <limits.h>


namespace mummer {
namespace mummer {

static const unsigned int BITADD[256] = { UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//0-9
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//10-19
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//20-29
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//30-39
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//40-49
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//50-59
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 0,        UINT_MAX, 1,        UINT_MAX, UINT_MAX,//60-69 65:A, 67:C
                                         UINT_MAX, 2,        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//70-79 71:G
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 3,        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//80-89 84:T
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 0,        UINT_MAX, 1,       //90-99 97:a, 99: c
                                         UINT_MAX, UINT_MAX, UINT_MAX, 2,        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//100-109 103:g
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 3,        UINT_MAX, UINT_MAX, UINT_MAX,//110-119 116:t
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//120-129
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//130-139
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//140-149
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//150-159
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//160-169
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//170-179
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//180-189
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//190-199
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//200-209
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//210-219
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//220-229
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//230-239
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//240-249
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX };//250-255

// Stores the LCP array in an unsigned char (0-255).  Values larger
// than or equal to 255 are stored in a sorted array.
// Simulates a vector<int> LCP;
struct vec_uchar {
  struct item_t{
    item_t(){}
    item_t(size_t i, int v) { idx = i; val = v; }
    size_t idx; int val;
    bool operator < (item_t t) const { return idx < t.idx;  }
  };
  std::vector<unsigned char> vec;  // LCP values from 0-65534
  std::vector<item_t> M;
  void resize(size_t N) { vec.resize(N); }
  // Vector X[i] notation to get LCP values.
  int operator[] (size_t idx) const {
    if(vec[idx] == std::numeric_limits<unsigned char>::max())
      return lower_bound(M.begin(), M.end(), item_t(idx,0))->val;
    else
      return vec[idx];
  }
  // Actually set LCP values, distingushes large and small LCP
  // values.
  void set(size_t idx, int v) {
    if(v >= std::numeric_limits<unsigned char>::max()) {
      vec.at(idx) = std::numeric_limits<unsigned char>::max();
      M.push_back(item_t(idx, v));
    }
    else { vec.at(idx) = (unsigned char)v; }
  }
  // Once all the values are set, call init. This will assure the
  // values >= 255 are sorted by index for fast retrieval.
  void init() { sort(M.begin(), M.end()); std::cerr << "M.size()=" << M.size() << std::endl; std::vector<item_t>(M).swap(M);}

  long index_size_in_bytes() const {
      long indexSize = 0L;
      indexSize += sizeof(vec) + vec.capacity()*sizeof(unsigned char);
      indexSize += sizeof(M) + M.capacity()*(sizeof(size_t)+sizeof(int));
      return indexSize;
  }
};

// Match find by findMEM.
struct match_t {
  match_t() { ref = 0; query = 0, len = 0; }
  match_t(long r, long q, long l) { ref = r; query = q; len = l; }
  long ref; // position in reference sequence
  long query; // position in query
  long len; // length of match
};

struct saTuple_t {
    saTuple_t(): left(0), right(0) {}
    saTuple_t(unsigned int l, unsigned int r): left(l), right(r){}
    unsigned int left;
    unsigned int right;
};

// depth : [start...end]
struct interval_t {
  interval_t() { start = 1; end = 0; depth = -1; }
  interval_t(long s, long e, long d) { start = s; end = e; depth = d; }
  void reset(long e) { start = 0; end = e; depth = 0; }
  long depth, start, end;
  long size() { return end - start + 1; }
};

struct sparseSA {
  const std::vector<std::string>& descr; // Descriptions of concatenated sequences.
  const std::vector<long>&        startpos; // Lengths of concatenated sequences.
  const long                      maxdescrlen; // Maximum length of the sequence description, used for formatting.
  const bool                      _4column; // Use 4 column output format.

  const long K;                 // suffix sampling, K = 1 every suffix, K = 2 every other suffix, K = 3, every 3rd sffix
  const std::string               S;  //!< Reference to sequence data.
  const long N;                       //!< Length of the sequence.
  const long                      logN; // ceil(log(N))
  const long                      NKm1; // N/K - 1
  std::vector<unsigned int> SA; // Suffix array.
  std::vector<int>          ISA; // Inverse suffix array.
  vec_uchar                 LCP; // Simulates a vector<int> LCP.
  std::vector<int>          CHILD; //child table
  std::vector<saTuple_t>    KMR;

  const bool hasChild;
  const bool hasSufLink;
  //fields for lookup table of sa intervals to a certain small depth
  const bool hasKmer;
  const long kMerSize;
  long       kMerTableSize;
  int        sparseMult;
  const bool printSubstring;
  const bool printRevCompForw;
  const bool nucleotidesOnly;

  long index_size_in_bytes(){
      long indexSize = 0L;
      indexSize += sizeof(printRevCompForw);
      indexSize += sizeof(printSubstring);
      indexSize += sizeof(sparseMult);
      indexSize += sizeof(hasSufLink);
      indexSize += sizeof(hasChild);
      indexSize += sizeof(K);
      indexSize += sizeof(NKm1);
      indexSize += sizeof(logN);
      indexSize += sizeof(N);
      indexSize += sizeof(_4column);
      indexSize += sizeof(maxdescrlen);
      indexSize += sizeof(descr);
      indexSize += sizeof(hasKmer);
      indexSize += sizeof(kMerSize);
      indexSize += sizeof(kMerTableSize);
      indexSize += sizeof(nucleotidesOnly);
      for(size_t i = 0; i < descr.size(); i++){
          indexSize += descr[i].capacity();
      }
      indexSize += sizeof(startpos) + startpos.capacity()*sizeof(long);
      indexSize += S.capacity();
      indexSize += sizeof(SA) + SA.capacity()*sizeof(unsigned int);
      indexSize += sizeof(ISA) + ISA.capacity()*sizeof(int);
      indexSize += sizeof(CHILD) + CHILD.capacity()*sizeof(int);
      indexSize += sizeof(KMR) + KMR.capacity()*(2*sizeof(unsigned int));
      indexSize += LCP.index_size_in_bytes();
      return indexSize;
  }

  // Maps a hit in the concatenated sequence set to a position in that sequence.
  void from_set(long hit, long &seq, long &seqpos) const {
    // Use binary search to locate index of sequence and position
    // within sequence.
    std::vector<long>::const_iterator it = upper_bound(startpos.begin(), startpos.end(), hit);   // SG: should use vector<long>::const_iterator
    seq = distance(startpos.begin(), it) - 1;
    it--;
    seqpos = hit - *it;
  }

  // Constructor builds sparse suffix array.
  sparseSA(const std::string& S_, const std::vector<std::string>& descr_, const std::vector<long>& startpos_, bool __4column,
           long K_, bool suflink_, bool child_, bool kmer_, int sparseMult_, int kMerSize_, bool printSubstring_,
           bool printRevCompForw_, bool nucleotidesOnly_);

  // Modified Kasai et all for LCP computation.
  void computeLCP();
  //Modified Abouelhoda et all for CHILD Computation.
  void computeChild();
  //build look-up table for sa intervals of kmers up to some depth
  void computeKmer();

  // Radix sort required to construct transformed text for sparse SA construction.
  void radixStep(int *t_new, int *SA, long &bucketNr, long *BucketBegin, long l, long r, long h);

  // Prints match to cout.
  void print_match(match_t m) const;
  void print_match(match_t m, std::vector<match_t> &buf) const; // buffered version
  void print_match(std::string meta, std::vector<match_t> &buf, bool rc) const; // buffered version

  // Binary search for left boundry of interval.
  inline long bsearch_left(char c, long i, long s, long e) const;
  // Binary search for right boundry of interval.
  inline long bsearch_right(char c, long i, long s, long e) const;

  // Simple suffix array search.
  inline bool search(std::string &P, long &start, long &end);

  // Simple top down traversal of a suffix array.
  inline bool top_down(char c, long i, long &start, long &end) const;
  inline bool top_down_faster(char c, long i, long &start, long &end) const;
  inline bool top_down_child(char c, interval_t &cur) const;

  // Traverse pattern P starting from a given prefix and interval
  // until mismatch or min_len characters reached.
  inline void traverse(std::string &P, long prefix, interval_t &cur, int min_len) const;
  inline void traverse_faster(const std::string &P,const long prefix, interval_t &cur, int min_len) const;

  // Simulate a suffix link.
  inline bool suffixlink(interval_t &m) const;

  // Expand ISA/LCP interval. Used to simulate suffix links.
  inline bool expand_link(interval_t &link) const {
    long thresh = 2 * link.depth * logN, exp = 0; // Threshold link expansion.
    long start = link.start;
    long end = link.end;
    while(LCP[start] >= link.depth) {
      exp++;
      if(exp >= thresh) return false;
      start--;
    }
    while(end < NKm1 && LCP[end+1] >= link.depth) {
      exp++;
      if(exp >= thresh) return false;
      end++;
    }
    link.start = start; link.end = end;
    return true;
  }

  // Given a position i in S, finds a left maximal match of minimum
  // length within K steps.
  inline void find_Lmaximal(std::string &P, long prefix, long i, long len, std::vector<match_t> &matches, int min_len, bool _forward, bool print) const;

  // Given an interval where the given prefix is matched up to a
  // mismatch, find all MEMs up to a minimum match depth.
  void collectMEMs(std::string &P, long prefix, interval_t mli, interval_t xmi, std::vector<match_t> &matches, int min_len, bool forward_, bool print) const;

  // Find all MEMs given a prefix pattern offset k.
  void findMEM(long k, std::string &P, std::vector<match_t> &matches, int min_len, bool forward_, bool print) const;

  // NOTE: min_len must be > 1
  void findMAM(std::string &P, std::vector<match_t> &matches, int min_len, long& memCount, bool forward_, bool print) const;
  // Returns true if the position p1 in the query pattern and p2 in
  // the reference is left maximal.
  inline bool is_leftmaximal(std::string &P, long p1, long p2) const {
    return p1 == 0 || p2 == 0 || P[p1 - 1] != S[p2 - 1];
  }

  // Maximal Almost-Unique Match (MAM). Match is unique in the indexed
  // sequence S. as computed by MUMmer version 2 by Salzberg
  // et. al. Note this is a "one-sided" query. It "streams" the query
  // P throught he index.  Consequently, repeats can occur in the
  // pattern P.
  void MAM(std::string &P, std::vector<match_t> &matches, int min_len, long& memCount, bool forward_, bool print) const {
    if(K != 1) return;  // Only valid for full suffix array.
    findMAM(P, matches, min_len, memCount, forward_, print);
  }

  // Find Maximal Exact Matches (MEMs)
  void MEM(std::string &P, std::vector<match_t> &matches, int min_len, bool print, long& memCount, bool forward_, int num_threads = 1)const ;

  // Maximal Unique Match (MUM)
  void MUM(std::string &P, std::vector<match_t> &unique, int min_len, long& memCount, bool forward_, bool print) const;

  //save index to files
  void save(const std::string &prefix);

  //load index from file
  bool load(const std::string &prefix);

  //construct
  void construct();
};

} // namespace mummer
} // namespace mummer

#endif // __sparseSA_hpp__

