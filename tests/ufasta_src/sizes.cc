#include <cstdlib>
#include <string>
#include <fstream>

#include "common.hpp"
#include <sizes_cmdline.hpp>

int sizes_main(int argc, char *argv[]) {
  sizes_cmdline args(argc, argv);
  if(args.file_arg.empty()) {
    args.file_arg.push_back("/dev/stdin");
    if(isatty(0))
      std::cerr << "Warning: reading from terminal" << std::endl;
  }

  std::string line;
  for(auto file : args.file_arg) {
    try {
      std::ifstream is;
      is.exceptions(std::ios::failbit|std::ios::badbit);
      is.open(file);
      if((args.file_arg.size() > 1 && !args.quiet_flag) || args.verbose_flag)
        std::cout << "==> " << file << " <==\n";

      int    c;
      size_t offset = 0;

      // Skip up to first header
      for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
        skip_line(is);
        offset += is.gcount();
      }

      while(c != EOF) {
        if(!args.header_flag) {
          skip_line(is);
          offset += is.gcount();
        } else {
          std::getline(is, line);
          offset    += line.size() + 1;
          auto last  = line.find_first_of(" \t");
          std::cout.write(line.c_str() + 1, (last == std::string::npos ? line.size() : last) - 1);
          std::cout << ' ';
        }
        const size_t save_offset = offset;

        size_t size = 0;
        for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
          skip_line(is);
          offset += is.gcount();
          size   += is.gcount() - 1;
        }
        std::cout << size;
        if(args.index_flag) std::cout << ' ' << save_offset << ' ' << offset;
        std::cout << '\n';
      }
    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
