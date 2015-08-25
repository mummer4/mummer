#ifndef __sparseSA_hpp__
#define __sparseSA_hpp__

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <limits.h>


namespace mummer {
namespace mummer {

static const unsigned int BITADD[256] = {
  UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//0-9
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
    item_t() = default;
    item_t(size_t i, int v) : idx(i), val(v) { }
    size_t idx; int val;
    bool operator < (item_t t) const { return idx < t.idx;  }
  };
  std::vector<unsigned char> vec; // LCP values from 0-65534
  std::vector<item_t>        M;
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
  void init() { sort(M.begin(), M.end()); } // std::vector<item_t>(M).swap(M);}

  long index_size_in_bytes() const {
      long indexSize = 0L;
      indexSize += sizeof(vec) + vec.capacity()*sizeof(unsigned char);
      indexSize += sizeof(M) + M.capacity()*(sizeof(size_t)+sizeof(int));
      return indexSize;
  }
};

// Match find by findMEM.
struct match_t {
  match_t() : ref(0), query(0), len(0) { }
  match_t(long r, long q, long l) : ref(r), query(q), len(l) { }
  long ref; // position in reference sequence
  long query; // position in query
  long len; // length of match
};


struct saTuple_t {
    saTuple_t(): left(0), right(0) {}
    saTuple_t(unsigned int l, unsigned int r): left(l), right(r) {}
    unsigned int left;
    unsigned int right;
};

// depth : [start...end]
struct interval_t {
  interval_t() : depth(-1), start(1), end(0) { }
  interval_t(long s, long e, long d) : depth(d), start(s), end(e) { }
  void reset(long e) { start = 0; end = e; depth = 0; }
  long depth, start, end;
  long size() const { return end - start + 1; }
};

struct bounded_string {
  const char* const s_;
  const size_t      al_; // actual length
  const size_t      l_;  // length rounded to K

  bounded_string(const char* s, size_t l, long K) :
    s_(s),
    al_(l),
    l_(l + K + (l % K != 0 ? K - (l % K) : 0))
  { }
  bounded_string(const std::string s, long K) : bounded_string(s.c_str(), s.size(), K) { }

  char operator[](size_t i) const {
    if(__builtin_expect(i < al_, 1))
      return s_[i];
    return '$';
  }
  size_t size() const noexcept { return l_; }
  size_t length() const noexcept { return l_; }
  size_t capacity() const noexcept { return l_; }
  std::string substr (size_t pos = 0, size_t len = std::string::npos) const {
    pos = std::min(pos, al_);
    return std::string(std::min(pos, al_), std::min(len, al_ - pos));
  }
};

struct sparseSA {
  const std::vector<std::string>& descr; // Descriptions of concatenated sequences.
  const std::vector<long>&        startpos; // Lengths of concatenated sequences.
  const long                      maxdescrlen; // Maximum length of the sequence description, used for formatting.
  const bool                      _4column; // Use 4 column output format.

  const long                K;  // suffix sampling, K = 1 every suffix, K = 2 every other suffix, K = 3, every 3rd sffix
  const bounded_string      S;  //!< Reference to sequence data.
  const long                N;  //!< Length of the sequence.
  const long                logN; // ceil(log(N))
  const long                NKm1; // N/K - 1
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

  static sparseSA create_auto(const std::string& S, const std::vector<std::string>& descr_, const std::vector<long>& startpos_,
                              int min_len, bool nucleotidesOnly_, int K = 1);

  // Modified Kasai et all for LCP computation.
  void computeLCP();
  //Modified Abouelhoda et all for CHILD Computation.
  void computeChild();
  //build look-up table for sa intervals of kmers up to some depth
  void computeKmer();

  // Radix sort required to construct transformed text for sparse SA construction.
  void radixStep(int *t_new, int *SA, long &bucketNr, long *BucketBegin, long l, long r, long h);

  // Prints match to cout.
  void print_match(std::ostream& os, match_t m) const;
  //  void print_match(match_t m, std::vector<match_t> &buf) const; // buffered version
  void print_match(std::ostream& os, std::string meta, bool rc) const; // buffered version

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
  void traverse(const std::string &P, long prefix, interval_t &cur, int min_len) const;
  void traverse_faster(const std::string &P,const long prefix, interval_t &cur, int min_len) const;

  // Simulate a suffix link.
  bool suffixlink(interval_t &m) const;

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
  template<typename Output>
  void find_Lmaximal(const std::string &P, long prefix, long i, long len, int min_len, bool _forward, Output out) const;

  // NOTE: min_len must be > 1
  template<typename Output>
  void findMAM_each(const std::string &P, int min_len, bool forward_, Output out) const;

  void findMAM(const std::string &P, int min_len, bool forward_, std::ostream& os) const {
    findMAM_each(P, min_len, forward_, [&](const match_t& m) { print_match(os, m); });
  }

  void findMAM(const std::string &P, int min_len, bool forward_, std::vector<match_t>& matches) const {
    findMAM_each(P, min_len, forward_, [&](const match_t& m) { matches.push_back(m); });
  }

  // Returns true if the position p1 in the query pattern and p2 in
  // the reference is left maximal.
  inline bool is_leftmaximal(const std::string &P, long p1, long p2) const {
    return p1 == 0 || p2 == 0 || P[p1 - 1] != S[p2 - 1];
  }

  // Maximal Almost-Unique Match (MAM). Match is unique in the indexed
  // sequence S. as computed by MUMmer version 2 by Salzberg
  // et. al. Note this is a "one-sided" query. It "streams" the query
  // P throught he index.  Consequently, repeats can occur in the
  // pattern P.
  void MAM(const std::string& P, int min_len, bool forward_, std::ostream& os) const {
    if(K != 1) return;  // Only valid for full suffix array.
    findMAM(P, min_len, forward_, os);
  }
  void MAM(const std::string& P, int min_len, bool forward_, std::vector<match_t>& matches) const {
    if(K != 1) return;  // Only valid for full suffix array.
    findMAM(P, min_len, forward_, matches);
  }


  // Given an interval where the given prefix is matched up to a
  // mismatch, find all MEMs up to a minimum match depth.
  template<typename Output>
  void collectMEMs_each(const std::string &P, long prefix, interval_t mli, interval_t xmi, int min_len, bool forward_,
                        Output out) const;

  // Find all MEMs given a prefix pattern offset k.
  template<typename Output>
  void findMEM_k_each(const std::string &P, long k, int min_len, bool forward_, Output out) const;

  // Find Maximal Exact Matches (MEMs)
  template<typename Output>
  void findMEM_each(const std::string &P, int min_len, bool forward_, Output out) const {
    for(int k = 0; k < K; k++)
      findMEM_k_each(P, k, min_len, forward_, out);
  }

  void MEM(const std::string &P, int min_len, bool forward_, std::ostream& os) const {
    findMEM_each(P, min_len, forward_, [&](const match_t& m) { print_match(os, m); });
  }

  void MEM(const std::string &P, int min_len, bool forward_, std::vector<match_t>& matches) const {
    findMEM_each(P, min_len, forward_, [&](const match_t& m) { matches.push_back(m); });
  }

  // Maximal Unique Match (MUM)
  template<typename Output>
  void findMUM_each(const std::string &P, int min_len, bool forward_, Output out) const;

  void MUM(const std::string &P, int min_len, bool forward_, std::ostream& os) const {
    findMUM_each(P, min_len, forward_, [&](const match_t& m) { print_match(os, m); });
  }

  void MUM(const std::string &P, int min_len, bool forward_, std::vector<match_t>& matches) const {
    findMUM_each(P, min_len, forward_, [&](const match_t& m) { matches.push_back(m); });
  }

  //save index to files
  void save(const std::string &prefix);

  //load index from file
  bool load(const std::string &prefix);

  //construct
  void construct();
};


//
// Implementation of templated methods
//

// Finds maximal almost-unique matches (MAMs) These can repeat in the
// given query pattern P, but occur uniquely in the indexed reference S.
template<typename Output>
void sparseSA::findMAM_each(const std::string &P, int min_len, bool forward_, Output out) const {
  const long Plength = P.length();
  interval_t cur(0, N-1, 0);
  long       prefix  = 0;

  while(prefix < Plength) {
    // Traverse SA top down until mismatch or full string is matched.
    if(hasChild)
        traverse_faster(P, prefix, cur, Plength);
    else
        traverse(P, prefix, cur, Plength);
    if(cur.depth <= 1) { cur.depth = 0; cur.start = 0; cur.end = N-1; prefix++; continue; }
    if(cur.size() == 1 && cur.depth >= min_len) {
      if(is_leftmaximal(P, prefix, SA[cur.start])) {
	// Yes, it's a MAM.
	match_t m; m.ref = SA[cur.start]; m.query = prefix; m.len = cur.depth;
        if(printRevCompForw && !forward_) m.query = Plength-1-prefix;
        out(m);
      }
    }
    do {
      cur.depth = cur.depth-1;
      cur.start = ISA[SA[cur.start] + 1];
      cur.end   = ISA[SA[cur.end] + 1];
      prefix++;
      if( cur.depth == 0 || !expand_link(cur) ) { cur.depth = 0; cur.start = 0; cur.end = N-1; break; }
    } while(cur.depth > 0 && cur.size() == 1);
  }
  //  currentCount = memCount;
}

// Maximal Unique Match (MUM)
template<typename Output>
void sparseSA::findMUM_each(const std::string &P, int min_len, bool forward_, Output out) const {
  // Find unique MEMs.
  std::vector<match_t> matches;
  MAM(P, min_len, forward_, matches);
  //  memCount=0;

  struct by_ref {
    bool operator() (const match_t &a, const match_t &b) const {
      return (a.ref == b.ref) ? a.len > b.len : a.ref < b.ref;
    }
  };

  // Adapted from Stephan Kurtz's code in cleanMUMcand.c in MUMMer v3.20.
  long currentright, dbright         = 0;
  bool ignorecurrent, ignoreprevious = false;
  sort(matches.begin(), matches.end(), by_ref());
  for(long i = 0; i < (long)matches.size(); i++) {
    ignorecurrent = false;
    currentright = matches[i].ref + matches[i].len - 1;
    if(dbright > currentright)
      ignorecurrent = true;
    else {
      if(dbright == currentright) {
	ignorecurrent = true;
	if(!ignoreprevious && matches[i-1].ref == matches[i].ref)
	  ignoreprevious = true;
      }
      else {
	dbright = currentright;
      }
    }
    if(i > 0 && !ignoreprevious)
      out(matches[i-1]);
    ignoreprevious = ignorecurrent;
  }
  if(!ignoreprevious && !matches.empty())
    out(matches.back());
}

// For a given offset in the prefix k, find all MEMs.
template<typename Output>
void sparseSA::findMEM_k_each(const std::string &P, long k, int min_len, bool forward_, Output out) const {
  if(k < 0 || k >= K) { std::cerr << "Invalid k " << k << " [0, " << K << "]" << std::endl; return; }
  // Offset all intervals at different start points.
  long prefix = k;
  interval_t mli(0,N/K-1,0); // min length interval
  interval_t xmi(0,N/K-1,0); // max match interval

  // Right-most match used to terminate search.
  const int min_lenK = min_len - (sparseMult*K-1);

  while( prefix <= (long)P.length() - min_lenK) {//BUGFIX: used to be "prefix <= (long)P.length() - (K-k0)"
    if(hasChild)
      traverse_faster(P, prefix, mli, min_lenK);    // Traverse until minimum length matched.
    else
      traverse(P, prefix, mli, min_lenK);    // Traverse until minimum length matched.

    if(mli.depth > xmi.depth) xmi = mli;
    if(mli.depth <= 1) { mli.reset(N/K-1); xmi.reset(N/K-1); prefix+=sparseMult*K; continue; }

    if(mli.depth >= min_lenK) {
      if(hasChild)
        traverse_faster(P, prefix, xmi, P.length()); // Traverse until mismatch.
      else
        traverse(P, prefix, xmi, P.length()); // Traverse until mismatch.
      collectMEMs_each(P, prefix, mli, xmi, min_len, forward_, out); // Using LCP info to find MEM length.
      // When using ISA/LCP trick, depth = depth - K. prefix += K.
      prefix+=sparseMult*K;
      if( !hasSufLink ){ mli.reset(N/K-1); xmi.reset(N/K-1); continue; }
      else{
        int i = 0;
        bool succes  = true;
        while(i < sparseMult && (succes = suffixlink(mli))){
          suffixlink(xmi);
          i++;
        }
        if(!succes){
          mli.reset(N/K-1); xmi.reset(N/K-1); continue;
        }
      }
    }
    else {
      // When using ISA/LCP trick, depth = depth - K. prefix += K.
      prefix+=sparseMult*K;
      if( !hasSufLink) { mli.reset(N/K-1); xmi.reset(N/K-1); continue; }
      else{
        int i = 0;
        bool succes  = true;
        while(i < sparseMult && (succes = suffixlink(mli))){
          i++;
        }
        if(!succes){
          mli.reset(N/K-1); xmi.reset(N/K-1); continue;
        }
        xmi = mli;
      }
    }
  }
}

// Use LCP information to locate right maximal matches. Test each for
// left maximality.
template<typename Output>
void sparseSA::collectMEMs_each(const std::string &P, long prefix, interval_t mli, interval_t xmi, int min_len, bool forward_,
                           Output out) const {
  // All of the suffixes in xmi's interval are right maximal.
  for(long i = xmi.start; i <= xmi.end; i++) find_Lmaximal(P, prefix, SA[i], xmi.depth, min_len, forward_, out);

  if(mli.start == xmi.start && mli.end == xmi.end) return;


  while(xmi.depth >= mli.depth) {
    // Attempt to "unmatch" xmi using LCP information.
    if(xmi.end+1 < N/K) xmi.depth = std::max(LCP[xmi.start], LCP[xmi.end+1]);
    else xmi.depth = LCP[xmi.start];

    // If unmatched XMI is > matched depth from mli, then examine rmems.
    if(xmi.depth >= mli.depth) {
      // Scan RMEMs to the left, check their left maximality..
      while(LCP[xmi.start] >= xmi.depth) {
	xmi.start--;
	find_Lmaximal(P, prefix, SA[xmi.start], xmi.depth, min_len, forward_, out);
      }
      // Find RMEMs to the right, check their left maximality.
      while(xmi.end+1 < N/K && LCP[xmi.end+1] >= xmi.depth) {
	xmi.end++;
	find_Lmaximal(P, prefix, SA[xmi.end], xmi.depth, min_len, forward_, out);
      }
    }
  }
}

// Finds left maximal matches given a right maximal match at position i.
template<typename Output>
void sparseSA::find_Lmaximal(const std::string &P, long prefix, long i, long len, int min_len, bool forward_, Output out) const {
  const long Plength = P.length();
  // Advance to the left up to K steps.
  for(long k = 0; k < sparseMult*K; k++) {
    // If we reach the end or a mismatch, and the match is long enough, print.
    if(prefix == 0 || i == 0 || P[prefix-1] != S[i-1]) {
      if(len >= min_len)
        out(match_t(i, (!printRevCompForw || forward_) ? prefix : Plength-1-prefix, len));
      return; // Reached mismatch, done.
    }
    prefix--; i--; len++; // Continue matching.
  }
}

} // namespace mummer
} // namespace mummer

#endif // __sparseSA_hpp__

