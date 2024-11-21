#include <cstdlib>
#include <string>
#include <fstream>

#include <common.hpp>
#include <head_cmdline.hpp>

// With a negative value, display all but the last N entries of a
// file. Store the entries in a circular buffer of size N. Any
// overflow get printed out.
static int head_negative(const head_cmdline& args) {
  int                res     = EXIT_SUCCESS;
  const int          entries = std::abs(args.entries_arg);
  std::vector<entry> cache(entries);

  for(const auto& file : args.file_arg) {
    try {
      std::ifstream is;
      is.exceptions(std::ios::failbit|std::ios::badbit);
      is.open(file);
      if((args.file_arg.size() > 1 && !args.quiet_flag) || args.verbose_flag)
        std::cout << "==> " << file << " <==\n";

      int c;
      // Display unchanged up to first header
      std::string line;
      for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
        std::getline(is, line);
        std::cout << line << '\n';
      }

      int cz = 0; // cache size
      for(int i = 0; c != EOF; cz = std::min(cz + 1, entries), i = (i + 1) % entries) {
        entry& ce = cache[i];
        if(cz == entries) { // Print with a delay (overflow)
          for(int j = 0; j < ce.size; ++j)
            std::cout << ce.lines[j] << '\n';
        }
        ce.size = 0;
        ce.add_line(is); // Replace entry
        for(c = is.peek(); c != '>' && c != EOF; c = is.peek())
          ce.add_line(is);
      }

    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      res = EXIT_FAILURE;
    }
  }

  return res;
}

// With a positive value, display the first N entries
static int head_positive(const head_cmdline& args) {
  int res = EXIT_SUCCESS;
  const auto line_limit = args.bytes_given ? std::numeric_limits<std::streamoff>::max() : (std::streamoff)args.entries_arg;
  const auto byte_limit = args.bytes_given ? (std::streamoff)args.bytes_arg : std::numeric_limits<std::streamoff>::max();

  std::string line;
  for(const auto& file : args.file_arg) {
    try {
      std::ifstream is;
      is.exceptions(std::ios::failbit|std::ios::badbit);
      is.open(file);
      if((args.file_arg.size() > 1 && !args.quiet_flag) || args.verbose_flag)
        std::cout << "==> " << file << " <==\n";
      int c;
      // Display unchanged up to first header
      for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
        std::getline(is, line);
        std::cout << line << '\n';
      }
      for(std::streamoff i = 0; c != EOF && i < line_limit && is.tellg() < byte_limit; ++i) {
        std::getline(is, line);
        std::cout << line << '\n';
        for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
          std::getline(is, line);
          std::cout << line << '\n';
        }
      }
    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      res = EXIT_FAILURE;
    }
  }

  return res;
}

int head_main(int argc, char *argv[]) {
  head_cmdline args(argc, argv);
  if(args.file_arg.empty()) {
    args.file_arg.push_back("/dev/stdin");
    if(isatty(0))
      std::cerr << "Warning: reading from terminal" << std::endl;
  }

  if(args.entries_arg == 0) return EXIT_SUCCESS;
  if(args.entries_arg > 0)
    return head_positive(args);
  else
    return head_negative(args);
}
