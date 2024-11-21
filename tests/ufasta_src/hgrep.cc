#include <string>
#include <fstream>
#include <limits>
#include <algorithm>
#include <boost/regex.hpp>

#include "common.hpp"
#include <hgrep_cmdline.hpp>


int hgrep_main(int argc, char *argv[]) {
  hgrep_cmdline args(argc, argv);
  if(args.file_arg.empty()) {
    args.file_arg.push_back("/dev/stdin");
    if(isatty(0))
      std::cerr << "Warning: reading from terminal" << std::endl;
  }

  boost::regex regexp;

  try {
    regexp = args.regexp_arg;
  } catch(boost::bad_expression& error) {
    std::cerr << "Invalid regular expression: " << error.what() << "\n"
              << "  " << args.regexp_arg << "\n"
              << std::string(1 + error.position(), ' ') << "^" << std::endl;
    return EXIT_FAILURE;
  }

  size_t max_count = args.max_count_given ? args.max_count_arg : std::numeric_limits<size_t>::max();
  std::string line;
  for(auto file : args.file_arg) {
    try {
      std::ifstream is;
      is.exceptions(std::ios::failbit|std::ios::badbit);
      is.open(file);
      bool          skip = true;
      int           c;
      // Display unchanged up to first header
      for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
        std::getline(is, line);
        std::cout << line << '\n';
      }

      while(c != EOF) {
        std::getline(is, line);
        skip = !regex_search(++line.cbegin(), line.cend(), regexp) ^ args.invert_match_flag;
        if(skip) {
          for(c = is.peek(); c != '>' && c != EOF; c = is.peek())
            skip_line(is);
        } else {
          if(max_count-- == 0) break;
          std::cout << line << '\n';
          for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
            std::getline(is, line);
            std::cout << line << '\n';
          }
        }
      }
    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
