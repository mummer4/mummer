#ifndef __COMMON_H__
#define __COMMON_H__

#include <string.h>

#include <istream>
#include <fstream>
#include <limits>
#include <vector>
#include <algorithm>
#include <functional>
#include <array>
#include <random>

inline std::istream& skip_line(std::istream& is) {
  return is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

struct entry {
  int                      size;
  std::vector<std::string> lines;
  void add_line(std::istream& is) {
    if(size >= (int)lines.size())
      lines.push_back("");
    std::getline(is, lines[size]);
    ++size;
  }
};

// Return a pointer to the nb-th space separated token in str. str is
// not modified. Returns NULL if less than nb columns.
const char* find_token(uint32_t nb, const char* str);

// Seed a random generator
template <typename EngineT, std::size_t StateSize = EngineT::state_size>
void
seed_prg(EngineT& engine, const char* save = nullptr, const char* load = nullptr)
{
  using          engine_type    = typename EngineT::result_type;
  using          device_type    = std::random_device::result_type;
  using          seedseq_type   = std::seed_seq::result_type;
  constexpr auto bytes_needed   = StateSize * sizeof(engine_type);
  constexpr auto numbers_needed = (sizeof(device_type) < sizeof(seedseq_type))
    ? (bytes_needed / sizeof(device_type))
    : (bytes_needed / sizeof(seedseq_type));
  std::array<device_type, numbers_needed> numbers {};

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
  std::seed_seq seedseq(numbers.cbegin(), numbers.cend());
  engine.seed(seedseq);

  if(save) {
    std::ofstream os(save);
    for(size_t i = 0; i < numbers_needed; ++i)
      os << numbers[i] << '\n';
    if(!os.good())
      throw std::runtime_error(std::string("Failed writing seed to '") + save + "'");
  }
}

template<typename EngineT>
EngineT seeded_prg(const char* save = nullptr, const char* load = nullptr)
{
  EngineT res;
  seed_prg(res, save, load);
  return res;
}

size_t append_line(std::istream& is, std::vector<char>& buf, size_t off);
#endif /* __COMMON_H__ */
