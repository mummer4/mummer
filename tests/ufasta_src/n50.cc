#include <cstdlib>
#include <cassert>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "common.hpp"
#include <n50_cmdline.hpp>

static std::vector<uint32_t> get_nsizes(const std::vector<uint32_t>& in_sizes, bool others) {
  std::vector<uint32_t> res;
  if(in_sizes.empty() && !others) {
    res.push_back(50);
  } else {
    for(auto x : in_sizes) {
      if(x == 0 || x >= 100)
        n50_cmdline::error() << "Invalid Nx size '" << x << "': x must be in (0, 100)";
      res.push_back(x);
    }
  }
  std::sort(res.begin(), res.end());

  return res;
}

void update_stats(size_t size, std::vector<size_t>& sizes, size_t& total_size, double& E) {
  sizes.push_back(size);
  size_t o_total_size  = total_size;
  total_size          += size;
  E = ((double)o_total_size / (double)total_size) * E + (double)(size * size) / total_size;
}

size_t read_from_fasta(std::istream& is, std::vector<size_t>& sizes, size_t& total_size, double& E) {
  size_t      contig_i = 0;
  int         c;
  std::string line;

  // Skip up to first header
  for(c = is.peek(); c != '>' && c != EOF; c = is.peek())
    skip_line(is);

  while(c != EOF) {
    skip_line(is); // Ignore header

    size_t size = 0;
    for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
      std::getline(is, line);
      size += line.size();
    }
    ++contig_i;
    update_stats(size, sizes, total_size, E);
  }
  return contig_i;
}

size_t read_from_sizes(std::istream& is, std::vector<size_t>& sizes, size_t& total_size, double& E) {
  size_t size;
  size_t contig_i = 0;
  std::string line;

  while(is.peek() != EOF) {
    std::getline(is, line);
    size = std::stoul(line);
    ++contig_i;
    update_stats(size, sizes, total_size, E);
  }
  return contig_i;
}

int n50_main(int argc, char *argv[]) {
  n50_cmdline args(argc, argv);
  if(args.file_arg.empty()) {
    args.file_arg.push_back("/dev/stdin");
    if(isatty(0))
      std::cerr << "Warning: reading from terminal" << std::endl;
  }


  const std::vector<uint32_t> nsizes = get_nsizes(args.N_arg, !args.all_flag && (args.Esize_flag || args.sum_flag));

  std::vector<size_t> sizes;
  size_t              total_size = 0;
  double              E          = 0.0;
  size_t              contig_i   = 0;
  for(auto file : args.file_arg) {
    try {
      std::ifstream is;
      is.exceptions(std::ios::failbit|std::ios::badbit);
      is.open(file);

      if(args.from_sizes_flag)
        contig_i += read_from_sizes(is, sizes, total_size, E);
      else
        contig_i += read_from_fasta(is, sizes, total_size, E);
    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      return EXIT_FAILURE;
    } catch(std::runtime_error& e) {
      std::cerr << "Error with file '" << file << "': " << e.what() << std::endl;
    }
  }

  // Compute statistics
  std::sort(sizes.begin(), sizes.end(), [](size_t x, size_t y) { return x > y; });
  const size_t sum_size = total_size;
  if(args.size_given)
    total_size = args.size_arg;

  std::vector<size_t> sizes_for_n;
  for(auto s : nsizes) sizes_for_n.push_back(std::ceil((double)s * total_size / 100.0));
  assert(sizes_for_n.size() == nsizes.size());

  size_t csize = 0;
  size_t nsize = 0;
  for(auto s : sizes) {
    csize += s;
    for( ; nsize < nsizes.size() && csize >= sizes_for_n[nsize]; ++nsize) {
      if(!args.no_header_flag)
        std::cout << 'N' << nsizes[nsize] << ' ';
      std::cout << s << '\n';
    }
  }
  for( ; nsize < nsizes.size(); ++nsize) {
    if(!args.no_header_flag)
      std::cout << 'N' << nsizes[nsize] << ' ';
    std::cout << "-\n";
  }
  if(args.sum_flag || args.all_flag) {
    if(!args.no_header_flag)
      std::cout << (args.terse_flag ? "S" :"Sequence") << ' ';
    std::cout << sum_size << '\n';
  }
  if(args.average_flag || args.all_flag) {
    if(!args.no_header_flag)
      std::cout << (args.terse_flag ? "A" : "Average") << ' ';
    std::cout << ((double)sum_size / sizes.size()) << "\n";
  }
  if(args.Esize_flag || args.all_flag) {
    if(!args.no_header_flag)
      std::cout << (args.terse_flag ? "E" : "E-size") << ' ';
    std::cout << E << '\n';
  }
  if(args.count_flag || args.all_flag) {
    if(!args.no_header_flag)
      std::cout << (args.terse_flag ? "C" : "Count") << ' ';
    std::cout << contig_i << '\n';
  }

  return EXIT_SUCCESS;
}
