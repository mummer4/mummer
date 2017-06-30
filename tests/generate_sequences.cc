#include <random>
#include <chrono>
#include <fstream>
#include <string>
#include <algorithm>
#include "generate_sequences_cmdline.hpp"

std::mt19937 rand_gen;

template<typename T>
T min(T x, T y) {
  return std::min(x, y);
}
template<typename T, typename... Ts>
T min(T x, T y, Ts... zs) {
  return min(std::min(x, y), zs...);
}

static char bases[4] = { 'a', 'c', 'g', 't' };
unsigned int sesab(char c) {
  switch(c) {
  case 'a': return 0;
  case 'c': return 1;
  case 'g': return 2;
  case 't': return 3;
  }
  return 0;
}
char comp(char c) {
  switch(c) {
  case 'a': return 't';
  case 'c': return 'g';
  case 'g': return 'c';
  case 't': return 'a';
  default:
    return 'n';
  }
}


std::string sequence(size_t len) {
  std::uniform_int_distribution<int> dist(0, 3);
  std::string                        result;
  for(size_t i = 0; i < len; ++i)
    result += bases[dist(rand_gen)];
  return result;
}

void write_sequence(std::ostream& os, const std::string& seq) {
  static const size_t len = 70;
  for(size_t i = 0; i < seq.size(); i += len)
    os << seq.substr(i, len) << '\n';
}

void rev_comp(std::string& s) {
  size_t i = 0, j = s.size() - 1;
  for( ; i < j; ++i, --j) {
    char t = comp(s[i]);
    s[i] = comp(s[j]);
    s[j] = t;
  }
  if(i == j)
    s[i] = comp(s[i]);
}

void generate_reads(const std::string& genome, size_t number, size_t length, unsigned int error, std::ostream& os) {
  std::string errors;
  std::uniform_int_distribution<unsigned int> has_error(0, 99);
  std::uniform_int_distribution<unsigned int> error_type(0, 2);
  std::uniform_int_distribution<size_t> position(0, genome.size() - length);
  std::uniform_int_distribution<unsigned int> subst(1, 3);
  std::uniform_int_distribution<unsigned int> insertion(0, 3);
  std::uniform_int_distribution<unsigned int> reverse(0, 1);

  for(size_t i = 0; i < number; ++i) {
    errors.clear();
    std::string seq = genome.substr(position(rand_gen), length);
    if(reverse(rand_gen))
      rev_comp(seq);
    for(size_t j = 0; j < seq.size(); ++j) {
      if(has_error(rand_gen) < error) {
        switch(error_type(rand_gen)) {
        case 0: // substitution
          seq[j] = bases[(sesab(seq[j]) + subst(rand_gen)) % 4];
          errors += "s" + std::to_string(j);
          break;
        case 1: // insertion
          seq.insert(j, 1, bases[insertion(rand_gen)]);
          errors += "i" + std::to_string(j);
          break;
        case 2: // deletion
          seq.erase(j, 1);
          errors += "d" + std::to_string(j);
          break;
        }
      }
    }
    os << '>' << i << ' ' << errors << '\n';
    write_sequence(os, seq);
  }
}

int main(int argc, char *argv[]) {
  cmdline args(argc, argv);

  if(args.seed_file_given) {
    std::ifstream is(args.seed_file_arg);
    if(!is.good())
      cmdline::error() << "Failed to open seed file '" << args.seed_file_arg << "'";
    is >> rand_gen;
  } else if(args.seed_given) {
    rand_gen.seed(args.seed_arg);
  } else {
    rand_gen.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  }

  { std::string name = args.prefix_arg + "_seed";
    std::ofstream os(name);
    if(!os.good())
      cmdline::error() << "Failed to open seed file '" << name << "'";
    os << rand_gen;
    if(!os.good())
      cmdline::error() << "Failed to write to seed file '" << name << "'";
  }

  if(!args.errors_given)
    args.errors_arg = { 1, 5 };
  if(!args.lengths_given)
    args.lengths_arg = { 100, 1000 };
  if(!args.numbers_given)
    args.numbers_arg = { 100, 100 };


  std::string genome = sequence(args.genome_size_arg);
  { std::string name = args.prefix_arg + "_genome.fa";
    std::ofstream os(name);
    if(!os.good())
      cmdline::error() << "Failed to open genome file '" << name << "'";
    os << ">genome\n";
    write_sequence(os, genome);
    if(!os.good())
      cmdline::error() << "Failed to write to genome file '" << name << "'";
  }

  for(size_t i = 0; i < min(args.errors_arg.size(), args.lengths_arg.size(), args.numbers_arg.size()); ++i) {
    std::string name = args.prefix_arg + "_reads_" + std::to_string(i) + ".fa";
    std::ofstream os(name);
    if(!os.good())
      cmdline::error() << "Failed to open reads fasta file '" << name << "'";
    generate_reads(genome, args.numbers_arg[i], args.lengths_arg[i], args.errors_arg[i], os);
    if(!os.good())
      cmdline::error() << "Failed to write to reads fasta file '" << name << "'";
  }

  return 0;
}
