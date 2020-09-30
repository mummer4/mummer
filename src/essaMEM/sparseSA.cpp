#include <math.h>
#include <pthread.h>
#include <limits.h>
#include <stdio.h>
#include <stack>
#include <assert.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include <mummer/sparseSA.hpp>
#include <mummer/sparseSA_imp.hpp>
#include <mummer/timer.hpp>
#include <compactsufsort/compactsufsort.hpp>

namespace mummer {
namespace mummer {

// LS suffix sorter (integer alphabet).
extern "C" { void suffixsort(int *x, int *p, int n, int k, int l); };

// Get maximum query sequence description length.
static size_t max_len(const std::vector<std::string>& descr) {
  size_t res = 0;
  for(size_t i = 0; i < descr.size(); ++i)
    res = std::max(res, descr[i].length());
  return res;
}

// static std::string prepare_sequence(std::string S, long K) {
//   // Increase string length so divisible by K.
//   // Don't forget to count $ termination character.
//   if(S.length() % K != 0) {
//     const long appendK = K - S.length() % K ;
//     for(long i = 0; i < appendK; i++) S += '$';
//   }
//   // Make sure last K-sampled characeter is this special character as well!!
//   for(long i = 0; i < K; i++) S += '$'; // Append "special" end character. Note: It must be lexicographically less.
//   return S;
// }

sparseSA::sparseSA(bounded_string&& S_,
                   bool __4column, long K_, bool suflink_, bool child_, bool kmer_, int sparseMult_,
                   int kMerSize_, bool nucleotidesOnly_)
  : sparseSA_aux(S_.length(), K_, __4column, child_, suflink_, kmer_, kMerSize_, sparseMult_, nucleotidesOnly_)
  , S(S_)
  , LCP(SA)
{ }

sparseSA sparseSA::create_auto(const char* S, size_t Slen, int min_len, bool nucleotidesOnly_, int K,
                               bool off48) {
  const bool suflink    = K < 4;
  const bool child      = K >= 4;
  int        sparseMult = 1;
  if(!suflink || child) {
    sparseMult = K >= 4
      ? (int) std::max((min_len-10)/K,1)
      : (int) std::max((min_len-12)/K,1);
  }
  const int kmer = std::max(0,std::min(10,min_len - sparseMult*K + 1));
  sparseSA res(S, Slen, true /* 4column */, K, suflink, child, kmer>0, sparseMult,
               kmer, nucleotidesOnly_);
  res.construct(off48);
  return res;
}

long sparseSA::index_size_in_bytes() const {
  throw std::runtime_error("TODO: broken");
      // long indexSize = 0L;
      // indexSize += sizeof(sparseMult);
      // indexSize += sizeof(hasSufLink);
      // indexSize += sizeof(hasChild);
      // indexSize += sizeof(K);
      // indexSize += sizeof(NKm1);
      // indexSize += sizeof(logN);
      // indexSize += sizeof(N);
      // indexSize += sizeof(_4column);
      // indexSize += sizeof(hasKmer);
      // indexSize += sizeof(kMerSize);
      // indexSize += sizeof(kMerTableSize);
      // indexSize += sizeof(nucleotidesOnly);
      // indexSize += S.capacity();
      // indexSize += sizeof(SA) + SA.capacity()*sizeof(unsigned int);
      // indexSize += sizeof(ISA) + ISA.capacity()*sizeof(int);
      // indexSize += sizeof(CHILD) + CHILD.capacity()*sizeof(int);
      // indexSize += sizeof(KMR) + KMR.capacity()*(2*sizeof(unsigned int));
      // indexSize += LCP.index_size_in_bytes();
      // return indexSize;
}

// Uses the algorithm of Kasai et al 2001 which was described in
// Manzini 2004 to compute the LCP array. Modified to handle sparse
// suffix arrays and inverse sparse suffix arrays.
void sparseSA::computeLCP() {
  TIME_FUNCTION;
  sparseSA_imp::computeLCP(LCP, S, SA, ISA, N, K);
}

// Child array construction algorithm
void sparseSA::computeChild() {
  //  TIME_FUNCTION;

  for(int i = 0; i < N/K; i++){
    CHILD[i] = -1;
  }
        //Compute up and down values
        int lastIndex  = -1;
        std::stack<int,std::vector<int> > stapelUD;
        stapelUD.push(0);
        for(int i = 1; i < N/K; i++){
            while(LCP[i] < LCP[stapelUD.top()]){
                lastIndex = stapelUD.top();
                stapelUD.pop();
                if(LCP[i] <= LCP[stapelUD.top()] && LCP[stapelUD.top()] != LCP[lastIndex]){
                    CHILD[stapelUD.top()] = lastIndex;
                }
            }
            //now LCP[i] >= LCP[top] holds
            if(lastIndex != -1){
                CHILD[i-1] = lastIndex;
                lastIndex = -1;
            }
            stapelUD.push(i);
        }
        while(0 < LCP[stapelUD.top()]){//last row (fix for last character of sequence not being unique
            lastIndex = stapelUD.top();
            stapelUD.pop();
            if(LCP[stapelUD.top()] != LCP[lastIndex]){
                CHILD[stapelUD.top()] = lastIndex;
            }
        }
        //Compute Next L-index values
        std::stack<int,std::vector<int> > stapelNL;
        stapelNL.push(0);
        for(int i = 1; i < N/K; i++){
            while(LCP[i] < LCP[stapelNL.top()])
                stapelNL.pop();
            lastIndex = stapelNL.top();
            if(LCP[i] == LCP[lastIndex]){
                stapelNL.pop();
                CHILD[lastIndex] = i;
            }
            stapelNL.push(i);
        }
}

// Look-up table construction algorithm
void sparseSA::computeKmer() {
  TIME_FUNCTION;

    std::stack<interval_t> intervalStack;
    std::stack<unsigned int> indexStack;

    interval_t curInterval(0,N/K-1,0);
    unsigned int curIndex = 0;
    unsigned int newIndex = 0;

    intervalStack.push(interval_t(0,N/K-1,0));
    indexStack.push(curIndex);

    while(!intervalStack.empty()){
        curInterval = intervalStack.top(); intervalStack.pop();
        curIndex = indexStack.top(); indexStack.pop();
        if(curInterval.depth == kMerSize){
            if(curIndex < kMerTableSize){
                KMR[curIndex].left = curInterval.start;
                KMR[curIndex].right = curInterval.end;
            }
        }
        else{
            if(hasChild){//similar to function traverse_faster
                //walk up to depth KMERSIZE or new child
                long curLCP; //max LCP of this interval
                if (curInterval.start == curInterval.end)
                    curLCP = N - SA[curInterval.start];
                else if (curInterval.start < CHILD[curInterval.end] && CHILD[curInterval.end] <= curInterval.end)
                    curLCP = LCP[CHILD[curInterval.end]];
                else
                    curLCP = LCP[CHILD[curInterval.start]];
                long minimum = std::min(curLCP,kMerSize);
                newIndex = curIndex;
                while(curInterval.depth < (long) minimum){
                    unsigned int character = S[SA[curInterval.start]+curInterval.depth];
                    newIndex = (newIndex << 2) | BITADD[character];
                    curInterval.depth ++;
                }
                if(curInterval.depth == kMerSize){//reached KMERSIZE in the middle of an edge
                    if(newIndex < kMerTableSize){
                        KMR[newIndex].left = curInterval.start;
                        KMR[newIndex].right = curInterval.end;
                    }
                }
                else{//find child intervals
                    long left = curInterval.start;
                    long right = CHILD[curInterval.end];
                    curIndex = newIndex;
                    if(curInterval.start >= right || right > curInterval.end)
                        right = CHILD[curInterval.start];
                    //now left and right point to first child
                    newIndex = (curIndex << 2) | BITADD[(int)S[SA[left]+curInterval.depth]];
                    if(newIndex < kMerTableSize){
                        intervalStack.push(interval_t(left,right-1,curInterval.depth+1));
                        indexStack.push(newIndex);
                    }
                    left = right;
                    //while has next L-index
                    while(CHILD[right] > right && LCP[right] == LCP[CHILD[right]]){
                        right = CHILD[right];
                        newIndex = (curIndex << 2) | BITADD[(int)S[SA[left]+curInterval.depth]];
                        if(newIndex < kMerTableSize){
                            intervalStack.push(interval_t(left,right-1,curInterval.depth+1));
                            indexStack.push(newIndex);
                        }
                        left = right;
                    }
                    //last interval
                    newIndex = (curIndex << 2) | BITADD[(int)S[SA[left]+curInterval.depth]];
                    if(newIndex < kMerTableSize){
                        intervalStack.push(interval_t(left,curInterval.end,curInterval.depth+1));
                        indexStack.push(newIndex);
                    }
                }
            }
            else{//similar to function traverse
                long start = curInterval.start; long end = curInterval.end;

                while(start <= curInterval.end){
                    unsigned int character = S[SA[start]+curInterval.depth];
                    newIndex = (curIndex << 2) | BITADD[character];
                    start = curInterval.start;
                    end = curInterval.end;
                    top_down_faster(character, curInterval.depth, start, end);
                    if(newIndex < kMerTableSize){
                        intervalStack.push(interval_t(start,end,curInterval.depth+1));
                        indexStack.push(newIndex);
                    }
                    // Advance to next interval.
                    start = end+1; end = curInterval.end;
                }
            }
        }
    }
}

bool vector_32_48::save(std::ostream&& os) const {
  size_t        size     = this->size();
  size_t        is_small = this->is_small;
  os.write((const char*)&size, sizeof(size));
  os.write((const char*)&is_small, sizeof(is_small));
  if(is_small) {
    os.write((const char*)small.data(), size * sizeof(int));
  } else {
    os.write((const char*)large.m_base32, size * sizeof(uint32_t));
    os.write((const char*)large.m_base16, size * sizeof(uint16_t));
  }
  return os.good();
}

bool vector_32_48::load(std::istream&& is) {
  size_t        size;
  size_t        is_small;
  is.read((char*)&size, sizeof(size));
  is.read((char*)&is_small, sizeof(is_small));
  resize(size, !is_small);
  if(is_small) {
    is.read((char*)small.data(), size * sizeof(int));
  } else {
    is.read((char*)large.m_base32, size * sizeof(uint32_t));
    is.read((char*)large.m_base16, size * sizeof(uint16_t));
  }
  return is.good();
}

bool sparseSA_aux::save(std::ostream&& os) const {
  os.write((const char*)&N,sizeof(N));
  os.write((const char*)&K,sizeof(K));
  os.write((const char*)&logN,sizeof(logN));
  os.write((const char*)&NKm1,sizeof(NKm1));
  os.write((const char*)&_4column,sizeof(_4column));
  os.write((const char*)&hasSufLink,sizeof(hasSufLink));
  os.write((const char*)&hasChild,sizeof(hasChild));
  os.write((const char*)&hasKmer,sizeof(hasKmer));
  os.write((const char*)&kMerSize,sizeof(kMerSize));
  os.write((const char*)&sparseMult,sizeof(sparseMult));
  os.write((const char*)&nucleotidesOnly,sizeof(nucleotidesOnly));
  return os.good();
}

bool sparseSA_aux::load(std::istream&& is) {
  is.read((char*)&N,sizeof(N));
  is.read((char*)&K,sizeof(K));
  is.read((char*)&logN,sizeof(logN));
  is.read((char*)&NKm1,sizeof(NKm1));
  is.read((char*)&_4column,sizeof(_4column));
  is.read((char*)&hasSufLink,sizeof(hasSufLink));
  is.read((char*)&hasChild,sizeof(hasChild));
  is.read((char*)&hasKmer,sizeof(hasKmer));
  is.read((char*)&kMerSize,sizeof(kMerSize));
  is.read((char*)&sparseMult,sizeof(sparseMult));
  is.read((char*)&nucleotidesOnly,sizeof(nucleotidesOnly));
  return is.good();
}

void vec_uchar::init() {
#ifdef _OPENMP
  openmp_qsort(M.begin(), M.end(), first_comp);
  assert(std::is_sorted(M.begin(), M.end(), first_comp));
#else
  sort(M.begin(), M.end(), first_comp);
#endif

  // Second, remove elements that are consecutive in a range
  size_t pidx = 0;
  large_type pval = 0;
  auto nend = std::remove_if(M.begin(), M.end(), [&](const item_t& a) -> bool {
      bool res = a.idx + a.val == pidx + pval && a.idx == pidx + 1;
      pidx = a.idx;
      pval = a.val;
      return res;
    });
  M.resize(nend - M.begin());
  M.shrink_to_fit();
  // Third, sort by beginning of compacted ranges
#ifdef _OPENMP
  openmp_qsort(M.begin(), M.end());
  assert(std::is_sorted(M.begin(), M.end()));
#else
  sort(M.begin(), M.end());
#endif
}

void vec_uchar::init_merge(const std::vector<item_vector>& Ms) {
  typedef item_vector::const_iterator iterator_type;
  struct elt {
    iterator_type it;
    iterator_type end;
  };
  auto comp = [&](const elt& a, const elt& b) -> bool {
    return first_comp(*b.it, *a.it);
  };
  std::vector<elt> heap;
  for(const auto& x : Ms)
    if(x.cbegin() != x.cend())
      heap.push_back(elt{ x.cbegin(), x.cend() });
  std::make_heap(heap.begin(), heap.end(), comp);

  size_t     pidx = 0;
  large_type pval = 0;
  while(!heap.empty()) {
    std::pop_heap(heap.begin(), heap.end(), comp);
    const auto  lend = heap.back().end;
    const auto& fit  = heap.front().it;
    auto&       lit  = heap.back().it;
    for( ; lit != lend && (first_comp(*lit, *fit) || *lit == *fit); ++lit) {
      if(lit->idx + lit->val != pidx + pval || lit->idx != pidx + 1)
        M.push_back(*lit);
      pidx = lit->idx;
      pval = lit->val;
    }
    if(lit != lend)
      std::push_heap(heap.begin(), heap.end(), comp);
    else
      heap.pop_back();
  }
  openmp_qsort(M.begin(), M.end());
}

bool vec_uchar::save(std::ostream&& os) const {
  const size_t sizeLCP = vec.size();
  const size_t sizeM   = M.size();
  os.write((const char*)&sizeLCP, sizeof(sizeLCP));
  os.write((const char*)&sizeM,   sizeof(sizeM));
  os.write((const char*)vec.data(), sizeLCP * sizeof(unsigned char));
  os.write((const char*)M.data(),   sizeM * sizeof(vec_uchar::item_t));
  return os.good();
}

bool vec_uchar::load(std::istream&& is) {
  size_t  sizeLCP;
  size_t  sizeM;
  is.read((char*)&sizeLCP, sizeof(sizeLCP));
  is.read((char*)&sizeM,   sizeof(sizeM));
  vec.resize(sizeLCP);
  M.resize(sizeM);
  is.read((char*)vec.data(), sizeLCP * sizeof(unsigned char));
  is.read((char*)M.data(),   sizeM * sizeof(vec_uchar::item_t));
  return is.good();
}


bool sparseSA::save(const std::string &prefix) const {
  //print auxiliary information
  if(!sparseSA_aux::save(prefix + ".aux"))
    return false;
  //print sa
  if(!SA.save(prefix + ".sa"))
    return false;
  if(!LCP.save(prefix + ".lcp"))
    return false;
  if(hasSufLink && !ISA.save(prefix + ".isa")) //print ISA if nec
    return false;
  if(hasChild){ //print child if nec
    const std::string child = prefix + ".child";
    std::ofstream child_s (child.c_str(), std::ios::binary);
    unsigned int sizeCHILD = CHILD.size();
    child_s.write((const char*)&sizeCHILD,sizeof(sizeCHILD));
    child_s.write((const char*)&CHILD[0],sizeCHILD*sizeof(int));
    if(!child_s.good()) return false;
  }
  if(hasKmer){ //print kmer if nec
    const std::string kmer = prefix + ".kmer";
    std::ofstream kmer_s (kmer.c_str(), std::ios::binary);
    unsigned int sizeKMR = KMR.size();
    kmer_s.write((const char*)&sizeKMR,sizeof(sizeKMR));
    kmer_s.write((const char*)&KMR[0],sizeKMR*sizeof(saTuple_t));
    if(!kmer_s.good()) return false;
  }
  return true;
}

bool sparseSA::load(const std::string &prefix) {
  // Load auxiliary infomation
  if(!sparseSA_aux::load(prefix + ".aux"))
    return false;
  //read sa
  if(!SA.load(prefix + ".sa"))
    return false;
  //read LCP
  LCP.sa = &SA;
  if(!LCP.load(prefix + ".lcp"))
    return false;
  if(hasSufLink && !ISA.load(prefix + ".isa")) //read ISA if nec
    return false;
  if(hasChild){ //read child if nec
    const std::string child = prefix + ".child";
    std::ifstream     child_s (child.c_str(), std::ios::binary);
    unsigned int      sizeCHILD;
    child_s.read((char*)&sizeCHILD,sizeof(sizeCHILD));
    CHILD.resize(sizeCHILD);
    child_s.read((char*)&CHILD[0],sizeCHILD*sizeof(int));
    if(!child_s.good()) return false;
  }
  if(hasKmer){ //read kmer table if nec
    const std::string kmer = prefix + ".kmer";
    std::ifstream     kmer_s (kmer.c_str(), std::ios::binary);
    unsigned int      sizeKMR;
    kmer_s.read((char*)&sizeKMR,sizeof(sizeKMR));
    KMR.resize(sizeKMR);
    kMerTableSize = sizeKMR;
    kmer_s.read((char*)&KMR[0],sizeKMR*sizeof(saTuple_t));
    if(!kmer_s.good()) return false;
  }

  return true;
}

void sparseSA::construct(bool off48){
  //  TIME_FUNCTION;

    if(K > 1) {
        // long bucketNr = 1;
        // int *intSA = new int[N/K+1];
        // for(int i = 0; i < N/K; i++) intSA[i] = i; // Init SA.
        // int* t_new = new int[N/K+1];
        // long* BucketBegin = new long[256]; // array to save current bucket beginnings
        // radixStep(t_new, intSA, bucketNr, BucketBegin, 0, N/K-1, 0); // start radix sort
        // t_new[N/K] = 0; // Terminate new integer string.
        // delete[] BucketBegin;

        throw "Not supported yet";
        // Suffix sort integer text.
        // std::cerr << "# suffixsort()" << std::endl;
        // suffixsort(t_new, intSA, N/K, bucketNr, 0);
        // std::cerr << "# DONE suffixsort()" << std::endl;

        // delete[] t_new;

        // // Translate suffix array.
        // SA.resize(N/K);
        // for (long i=0; i<N/K; i++) SA[i] = (unsigned int)intSA[i+1] * K;
        // delete[] intSA;

        // // Build ISA using sparse SA.
        // ISA.resize(N/K);
        // for(long i = 0; i < N/K; i++) { ISA[SA[i]/K] = i; }
    }
    else {
      SA.resize(N, off48);
      ISA.resize(N, off48);
      if(SA.is_small) {
        compactsufsort::create((const unsigned char*)(S + 0), (int*)SA.small.data(), N);
//#pragma omp parallel for
        for(long i = 0; i < N; ++i) { ISA.small[SA.small[i]] = i; }
      } else {
        compactsufsort::create((const unsigned char*)(S + 0), SA.large.begin(), N);
//#pragma omp parallel for
        for(long i = 0; i < N; ++i) { ISA.large[SA.large[i]] = i; }
      }
    }

    LCP.resize(N/K);
    // Use algorithm by Kasai et al to construct LCP array.
    computeLCP();  // SA + ISA -> LCP
    if(!hasSufLink){
      //ISA.clear(); // TODO: clear in vector32_48
    }
    if(hasChild){
        CHILD.resize(N/K);
        //Use algorithm by Abouelhoda et al to construct CHILD array
        computeChild();
    }
    if(hasKmer){
        kMerTableSize = 1 << (2*kMerSize);
        // std::cerr << "kmer table size: " << kMerTableSize << std::endl;
        KMR.resize(kMerTableSize, saTuple_t());
        //Use algorithm by Abouelhoda et al to construct CHILD array
        computeKmer();
    }

    //    NKm1 = N/K-1;

}

// Implements a variant of American flag sort (McIlroy radix sort).
// Recurse until big-K size prefixes are sorted. Adapted from the C++
// source code for the wordSA implementation from the following paper:
// Ferragina and Fischer. Suffix Arrays on Words. CPM 2007.
// void sparseSA::radixStep(int *t_new, int *SA, long &bucketNr, long *BucketBegin, long l, long r, long h) {
//   if(h >= K) return;
//   // first pass: count
//   std::vector<long> Sigma(256, 0); // Sigma counts occurring characters in bucket
//   for (long i = l; i <= r; i++) Sigma[ S[ SA[i]*K + h ] ]++; // count characters
//   BucketBegin[0] = l; for (long i = 1; i < 256; i++) { BucketBegin[i] = Sigma[i-1] + BucketBegin[i-1]; } // accumulate

//   // second pass: move (this variant does *not* need an additional array!)
//   unsigned char currentKey = 0;    // character of current bucket
//   long end = l-1+Sigma[currentKey]; // end of current bucket
//   long pos = l;                     // 'pos' is current position in bucket
//   while (1) {
//     if (pos > end) { // Reached the end of the bucket.
//       if (currentKey == 255) break; // Last character?
//       currentKey++; // Advance to next characer.
//       pos = BucketBegin[currentKey]; // Next bucket start.
//       end += Sigma[currentKey]; // Next bucket end.
//     }
//     else {
//       // American flag sort of McIlroy et al. 1993. BucketBegin keeps
//       // track of current position where to add to bucket set.
//       int tmp = SA[ BucketBegin[(int) S[ SA[pos]*K + h ] ] ];
//       SA[ BucketBegin[(int) S[ SA[pos]*K + h] ]++ ] = SA[pos];  // Move bucket beginning to the right, and replace
//       SA[ pos ] = tmp; // Save value at bucket beginning.
//       if (S[ SA[pos]*K + h ] == currentKey) pos++; // Advance to next position if the right character.
//     }
//   }
//   // recursively refine buckets and calculate new text:
//   long beg = l; end = l-1;
//   for (long i = 1; i < 256; i++) { // step through Sigma to find bucket borders
//     end += Sigma[i];
//     if (beg <= end) {
//       if(h == K-1) {
// 	for (long j = beg; j <= end; j++) {
// 	  t_new[ SA[j] ] = bucketNr; // set new text
// 	}
// 	bucketNr++;
//       }
//       else {
// 	radixStep(t_new, SA, bucketNr, BucketBegin, beg, end, h+1); // recursive refinement
//       }
//       beg = end + 1; // advance to next bucket
//     }
//   }
// }

// Binary search for left boundry of interval.
long sparseSA::bsearch_left(char c, long i, long s, long e) const {
  if(c == S[SA[s]+i]) return s;
  long l = s, r = e;
  while (r - l > 1) {
    long m = (l+r) / 2;
    if (c <= S[SA[m] + i]) r = m;
    else l = m;
  }
  return r;
}

// Binary search for right boundry of interval.
long sparseSA::bsearch_right(char c, long i, long s, long e) const {
  if(c == S[SA[e]+i]) return e;
  long l = s, r = e;
  while (r - l > 1) {
    long m = (l+r) / 2
;
    if (c < S[SA[m] + i]) r = m;
    else l = m;
  }
  return l;
}


// Simple top down traversal of a suffix array.
bool sparseSA::top_down(char c, long i, long &start, long &end) const {
  if(c < S[SA[start]+i]) return false;
  if(c > S[SA[end]+i]) return false;
  long l = bsearch_left(c, i, start, end);
  long l2 = bsearch_right(c, i, start, end);
  start = l; end = l2;
  return l <= l2;
}

// Top down traversal of the suffix array to match a pattern.  NOTE:
// NO childtab as in the enhanced suffix array (ESA).
bool sparseSA::search(const char* P, size_t Plen, long &start, long &end) const {
  start = 0; end = N/K - 1;
  long i = 0;
  while(i < (long)Plen) {
    if(top_down(P[i], i, start, end) == false) {
      return false;
    }
    i++;
  }
  return true;
}


// Traverse pattern P starting from a given prefix and interval
// until mismatch or min_len characters reached.
void sparseSA::traverse(const char* P, size_t Plen, long prefix, interval_t &cur, int min_len) const {
  if(hasKmer && cur.depth == 0 && min_len >= kMerSize){//free match first bases
    if((size_t)(prefix + kMerSize) > Plen) return;
    unsigned int index = 0;
    for(long i = 0; i < kMerSize; i++) {
      index = (index << 2 ) | BITADD[(int)P[prefix + i]];
    }
    if(index < kMerTableSize && KMR[index].right>0){
      cur.depth = kMerSize;
      cur.start = KMR[index].left;
      cur.end = KMR[index].right;
    }
    else if(index < kMerTableSize || nucleotidesOnly){
      return;//this results in no found seeds where the first KMERSIZE bases contain a non-ACGT character
    }
  }
  if(cur.depth >= min_len) return;
  while(prefix+cur.depth < (long)Plen) {
    long start = cur.start; long end = cur.end;
    // If we reach a mismatch, stop.
    if(top_down_faster(P[prefix+cur.depth], cur.depth, start, end) == false) return;

    // Advance to next interval.
    cur.depth += 1; cur.start = start; cur.end = end;

    // If we reach min_len, stop.
    if(cur.depth == min_len) return;
  }
}

// Traverse pattern P starting from a given prefix and interval
// until mismatch or min_len characters reached.
// Uses the child table for faster traversal
void sparseSA::traverse_faster(const char* P, size_t Plen, const long prefix, interval_t &cur, int min_len) const {
  if(hasKmer && cur.depth == 0 && min_len >= kMerSize){//free match first bases
    unsigned int index = 0;
    for(long i = 0; i < kMerSize; i++)
      index = (index << 2 ) | BITADD[(int)P[prefix + i]];
    if(index < kMerTableSize && KMR[index].right>0){
      cur.depth = kMerSize;
      cur.start = KMR[index].left;
      cur.end = KMR[index].right;
    }
    else if(index < kMerTableSize || nucleotidesOnly){
      return;//this results in no found seeds where the first KMERSIZE bases contain a non-ACGT character
    }
  }
  if(cur.depth >= min_len) return;
  long c = prefix + cur.depth;
  bool intervalFound = (size_t)c < Plen;
  int curLCP;//check if this is correct for root interval (unlikely case)
  if(cur.start < CHILD[cur.end] && CHILD[cur.end] <= cur.end)
    curLCP = LCP[CHILD[cur.end]];
  else
    curLCP = LCP[CHILD[cur.start]];
  if(intervalFound && cur.size() > 1 && curLCP == cur.depth)
    intervalFound = top_down_child(P[c], cur);
  else if(intervalFound)
    intervalFound = P[c] == S[SA[cur.start]+cur.depth];
  bool mismatchFound = false;
  while(intervalFound && !mismatchFound &&
        (size_t)c < Plen && cur.depth < min_len){
    c++;
    cur.depth++;
    if(cur.start != cur.end){
      int childLCP;
      //calculate LCP of child node, which is now cur. the LCP value
      //of the parent is currently c - prefix
      if(cur.start < CHILD[cur.end] && CHILD[cur.end] <= cur.end)
        childLCP = LCP[CHILD[cur.end]];
      else
        childLCP = LCP[CHILD[cur.start]];
      int minimum = std::min(childLCP,min_len);
      //match along branch
      while(!mismatchFound && (size_t)c < Plen && cur.depth < minimum){
        mismatchFound = S[SA[cur.start]+cur.depth] != P[c];
        c++;
        cur.depth += !mismatchFound;
      }
      intervalFound = (size_t)c < Plen && !mismatchFound &&
                                  cur.depth < min_len && top_down_child(P[c], cur);
    }
    else{
      while(!mismatchFound && (size_t)c < Plen && cur.depth < min_len){
        mismatchFound = (size_t)(SA[cur.start]+cur.depth) >= S.length() ||
                                          S[SA[cur.start]+cur.depth] != P[c];
        c++;
        cur.depth += !mismatchFound;
      }
    }
  }
}
//finds the child interval of cur that starts with character c
//updates left and right bounds of cur to child interval if found, or returns
//cur if not found (also returns true/false if found or not)
bool sparseSA::top_down_child(char c, interval_t &cur) const {
  long left = cur.start;
  long right = CHILD[cur.end];
  if(cur.start >= right || right > cur.end)
    right = CHILD[cur.start];
  //now left and right point to first child
  if(S[SA[cur.start]+cur.depth] == c){
    cur.end = right-1;
    return true;
  }
  left = right;
  //while has next L-index
  while(CHILD[right] > right && LCP[right] == LCP[CHILD[right]]){
    right = CHILD[right];
    if(S[SA[left]+cur.depth] == c){
      cur.start = left; cur.end = right - 1;
      return true;
    }
    left = right;
  }
  //last interval
  if(S[SA[left]+cur.depth] == c){
    cur.start = left;
    return true;
  }
  return false;
}

// Given SA interval apply binary search to match character c at
// position i in the search string. Adapted from the C++ source code
// for the wordSA implementation from the following paper: Ferragina
// and Fischer. Suffix Arrays on Words. CPM 2007.
bool sparseSA::top_down_faster(char c, long i, long &start, long &end) const {
  long l, r, m, r2=end, l2=start, vgl;
  bool found = false;
  long cmp_with_first = (long)c - (long)S[SA[start]+i];
  long cmp_with_last = (long)c - (long)S[SA[end]+i];
  if(cmp_with_first < 0) {
    l = start+1; l2 = start; // pattern doesn't occur!
  }
  else if(cmp_with_last > 0) {
    l = end+1; l2 = end;
    // pattern doesn't occur!
  }
  else {
    // search for left border:
    l = start; r = end;
    if (cmp_with_first == 0) {
      found = true; r2 = r;
    }
    else {
      while (r - l > 1) {
	m = (l+r) / 2;
	vgl = (long)c - (long)S[SA[m] + i];
	if (vgl <= 0) {
	  if (!found && vgl == 0) {
	    found = true;
	    l2 = m; r2 = r; // search interval for right border
	  }
	  r = m;
	}
	else l = m;
      }
      l = r;
    }
    // search for right border (in the range [l2:r2])
    if (!found) {
      l2 = l - 1; // pattern not found => right border to the left of 'l'
    }
    if (cmp_with_last == 0) {
      l2 = end; // right border is the end of the array
    }
    else {
      while (r2 - l2 > 1) {
	m = (l2 + r2) / 2;
	vgl = (long)c - (long)S[SA[m] + i];
	if (vgl < 0) r2 = m;
	else l2 = m;
      }
    }
  }
  start = l;
  end = l2;
  return l <= l2;
}


// Suffix link simulation using ISA/LCP heuristic.
bool sparseSA::suffixlink(interval_t &m) const {
  m.depth -= K;
  if( m.depth <= 0) return false;
  m.start = ISA[SA[m.start] / K + 1];
  m.end = ISA[SA[m.end] / K + 1];
  return expand_link(m);
}


///////////////////
// sparseSAMatch //
///////////////////
sparseSAMatch::sparseSAMatch(const std::string& S_, const std::vector<std::string>& descr_, const std::vector<long>& startpos_,
                             bool __4column, long K_, bool suflink_, bool child_, bool kmer_,
                             int sparseMult_, int kMerSize_, bool printSubstring_,
                             bool nucleotidesOnly_)
  : sparseSA(S_, __4column, K_, suflink_, child_, kmer_, sparseMult_, kMerSize_, nucleotidesOnly_)
  , descr(descr_)
  , startpos(startpos_)
  , maxdescrlen(max_len(descr_))
  , printSubstring(printSubstring_)
{ }

long sparseSAMatch::index_size_in_bytes() const {
  long indexSize = sparseSA::index_size_in_bytes();
  indexSize += sizeof(maxdescrlen);
  indexSize += sizeof(descr);
  for(size_t i = 0; i < descr.size(); i++){
    indexSize += descr[i].capacity();
  }
  indexSize += sizeof(startpos) + startpos.capacity()*sizeof(long);
  //  indexSize += sizeof(printRevCompForw);
  indexSize += sizeof(printSubstring);

  return indexSize;
}

void sparseSAMatch::from_set(long hit, long &seq, long &seqpos) const {
  // Use binary search to locate index of sequence and position
  // within sequence.
  std::vector<long>::const_iterator it = upper_bound(startpos.begin(), startpos.end(), hit);   // SG: should use vector<long>::const_iterator
  seq = distance(startpos.begin(), it) - 1;
  it--;
  seqpos = hit - *it;
}

// Print results in format used by MUMmer v3.  Prints results
// 1-indexed, instead of 0-indexed.
void sparseSAMatch::print_match(std::ostream& os, match_t m) const {
  if(_4column == false) {
    os << std::setw(8) << (m.ref + 1) << "  "
       << std::setw(8) << (m.query + 1) << "  "
       << std::setw(8) << m.len << '\n';
  }
  else {
    long refseq=0, refpos=0;
    from_set(m.ref, refseq, refpos); // from_set is slow!!!
    os << "  " << std::left << std::setw(maxdescrlen + 1) << descr[refseq] << std::right << ' '
       << std::setw(8) << (refpos + 1) << "  "
       << std::setw(8) << (m.query + 1) << "  "
       << std::setw(8) << m.len << '\n';
  }
  if(printSubstring){
    if(m.len > 53) os << S.substr(m.ref, 53) << " . . .\n";
    else os << m.ref << '\n';
  }
}

void sparseSAMatch::print_match(std::ostream& os, std::string meta, bool rc) const {
  if(!rc) os << "> " << meta << '\n';
  else os << "> " << meta << " Reverse\n";
}

struct thread_data {
  std::vector<long>  Kvalues;   // Values of K this thread should process.
  const sparseSA    *sa;        // Suffix array + aux informaton
  int                min_len;   // Minimum length of match.
  std::string       *P;         // Query string.
  bool               forward_;
};

} // namespace mummer
} // namespace mummer
