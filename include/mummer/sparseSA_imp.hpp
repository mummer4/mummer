#ifndef __SPARSESA_IMP_H__
#define __SPARSESA_IMP_H__

#ifdef _OPENMP
#include <omp.h>
#endif

// Implementation of some sparseSA functions
namespace mummer {
namespace sparseSA_imp {

#ifndef _OPENMP
template<typename Map, typename Seq, typename Vec>
void computeLCP(Map& LCP, const Seq& S, const Vec& SA, const Vec& ISA, const long N, const long K) {
  long h = 0;
  for(long i = 0; i < N / K; ++i) {
    const long m = ISA[i];
    if(m > 0) {
      const long bj  = SA[m-1];
      const long bi = i * K;
      while(bi + h < N && bj + h < N && S[bi + h] == S[bj + h])  ++h;
      LCP.set(m, h); //LCP[m] = h;
    } else {
      LCP.set(m, 0); // LCP[m]=0;
    }
    h = std::max(0L, h - K);
  }
  LCP.init();
}

#else // _OPENMP

template<typename Map, typename Seq, typename Vec>
void computeLCP(Map& LCP, const Seq& S, const Vec& SA, const Vec& ISA, const long N, const long K) {
  //  const long chunk = std::min((long)10000000, (long)(N / omp_get_num_threads()));
  std::vector<typename Map::item_vector> Ms;

#pragma omp parallel
  {
    long                      h = 0;
    typename Map::item_vector tM; // M array for this thread

#pragma omp for schedule(static)
    for(long i = 0; i < N / K; ++i) {
      const long m = ISA[i];
      if(m > 0) {
        const long bj  = SA[m-1];
        const long bi = i * K;
        while(bi + h < N && bj + h < N && S[bi + h] == S[bj + h])  ++h;
        LCP.set(m, h, tM); //LCP[m] = h;
      } else {
        LCP.set(m, 0); // LCP[m]=0;
      }
      h = std::max(0L, h - K);
    }

    std::sort(tM.begin(), tM.end(), Map::first_comp);
#pragma omp critical
    Ms.push_back(std::move(tM));
  }

  LCP.init_merge(Ms);
}

#endif // _OPENMP

} // namespace sparseSA_imp
} // namespace mummer

#endif /* __SPARSESA_IMP_H__ */
