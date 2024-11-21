#ifndef __SEEDED_PRG_H__
#define __SEEDED_PRG_H__

#include <fstream>
#include <random>
#include <algorithm>
#include <array>

template<typename EngineT, std::size_t StateSize = EngineT::state_size>
class LongSeed {

  using                 engine_type    = typename EngineT::result_type;
  using                 device_type    = std::random_device::result_type;
  using                 seedseq_type   = std::seed_seq::result_type;
  constexpr static auto bytes_needed   = StateSize * sizeof(engine_type);
  constexpr static auto numbers_needed = (sizeof(device_type) < sizeof(seedseq_type))
    ? (bytes_needed / sizeof(device_type))
    : (bytes_needed / sizeof(seedseq_type));

  std::array<device_type, numbers_needed> numbers;

public:
  typedef EngineT type;

  LongSeed(const char* save = nullptr, const char* load = nullptr) {
    if(load) {
      std::ifstream is(load);
      size_t i = 0;
      for( ; is && i < numbers_needed; ++i)
        is >> numbers[i];
      if(i != numbers_needed)
        std::runtime_error(std::string("Failed loading seed from '") + load + "'");
    } else {
      std::random_device rnddev {};
      std::generate(numbers.begin(), numbers.end(), std::ref(rnddev));
    }

    if(save) {
      std::ofstream os(save);
      for(size_t i = 0; i < numbers_needed; ++i)
        os << numbers[i] << '\n';
      if(!os.good())
        throw std::runtime_error(std::string("Failed writing seed to '") + save + "'");
    }
  }

  void seed(EngineT& engine) {
    std::seed_seq seedseq(numbers.cbegin(), numbers.cend());
    engine.seed(seedseq);
  }
  inline void operator()(EngineT& engine) { seed(engine); }

  EngineT seeded() {
    EngineT res;
    seed(res);
    return res;
  }
};

#endif /* __SEEDED_PRG_H__ */
