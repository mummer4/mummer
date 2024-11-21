#include <cstdlib>
#include <string>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <common.hpp>
#include <tail_cmdline.hpp>

// With a positive number, display the last N lines. Store lines in a
// circular buffer. Display content of circular buffer when reached
// the end of the input.
static int tail_positive(const tail_cmdline& args, int entries) {
  int                res = EXIT_SUCCESS;
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
      // Read into cache
      for(int i = 0; c != EOF; cz = std::min(cz + 1, entries), i = (i + 1) % entries) {
        entry& ce = cache[i];
        ce.size = 0;
        ce.add_line(is);
        for(c = is.peek(); c != '>' && c != EOF; c = is.peek())
          ce.add_line(is);
      }

      // Output cache
      for(int i = 0; i < cz; ++i) {
        const entry& ce = cache[i];
        for(int j = 0; j < ce.size; ++j)
          std::cout << ce.lines[j] << '\n';
      }
    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      res = EXIT_FAILURE;
    }
  }

  return res;
}

static int tail_negative(const tail_cmdline& args, const std::streamoff entries, const std::streamoff bytes) {
  int res = EXIT_SUCCESS;

  for(const auto& file : args.file_arg) {
    try {
      std::ifstream is;
      is.exceptions(std::ios::failbit|std::ios::badbit);
      is.open(file);
      if((args.file_arg.size() > 1 && !args.quiet_flag) || args.verbose_flag)
        std::cout << "==> " << file << " <==\n";
      int c = is.peek();
      for(std::streamoff i = 1; c != EOF && i < entries && is.tellg() < bytes; ++i) {
        if(c == '>') {
          skip_line(is);
          c = is.peek();
        }
        for( ; c != '>' && c != EOF; c = is.peek())
          skip_line(is);
      }
      if(c != EOF)
        std::cout << is.rdbuf();
    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      res = EXIT_FAILURE;
    }
  }

  return res;
}

template<typename T>
static long parse_nb(const T& s) {
  if(s.empty()) return 0;
  long res = 0;
  try {
    res= s.as_long(true);
  } catch(std::runtime_error& e) {
    tail_cmdline::error() << e.what();
  }
  return (s[0] == '+') ? -res : res;
}

int tail_main(int argc, char *argv[]) {
  tail_cmdline args(argc, argv);
  if(args.file_arg.empty()) {
    args.file_arg.push_back("/dev/stdin");
    if(isatty(0))
      std::cerr << "Warning: reading from terminal" << std::endl;
  }
  std::streamoff entries = 0;
  std::streamoff bytes   = 0;
  if(!args.bytes_given) {
    bytes = std::numeric_limits<std::streamoff>::max();
    try { entries = parse_nb(args.entries_arg);
    } catch(std::invalid_argument&) {
      tail_cmdline::error() << "Invalid number of entries '" << args.entries_arg << '\'';
    }
  } else {
    entries = std::numeric_limits<std::streamoff>::max();
    try { bytes = parse_nb(args.bytes_arg);
    } catch(std::invalid_argument&) {
      tail_cmdline::error() << "Invalid number of bytes '" << args.bytes_arg << '\'';
    }
    if(bytes > 0)
      tail_cmdline::error() << "Invalid bytes offset '" << args.bytes_arg << "'. Should be +K";
    bytes = -bytes;
  }

  if(entries == 0) return EXIT_SUCCESS;
  if(args.bytes_given)
    return tail_negative(args, entries, bytes);
  else if(entries > 0)
    return tail_positive(args, entries);
  else
    return tail_negative(args, std::abs(entries), bytes);

  return EXIT_SUCCESS;
}
