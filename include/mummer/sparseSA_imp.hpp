#ifndef __SPARSESA_IMP_H__
#define __SPARSESA_IMP_H__

// Implementation of some sparseSA functions
namespace mummer {
namespace sparseSA_imp {

template<typename Map, typename Seq, typename Vec>
void computeLCP(Map& LCP, const Seq& S, const Vec& SA, const Vec& ISA, const long N, const long K) {
#pragma omp parallel
  {
    long h = 0;
#pragma omp for schedule(static)
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
  }
}

} // namespace sparseSA_imp
} // namespace mummer

#endif /* __SPARSESA_IMP_H__ */
