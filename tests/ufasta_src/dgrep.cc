#include <string>
#include <fstream>
#include <limits>
#include <algorithm>
#include <boost/regex.hpp>

#include "common.hpp"
#include <dgrep_cmdline.hpp>


int dgrep_main(int argc, char *argv[]) {
  dgrep_cmdline args(argc, argv);
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
  size_t n_line = 0;
  std::string line, header;
  for(auto file : args.file_arg) {
    try {
      std::ifstream is;
      is.exceptions(std::ios::failbit|std::ios::badbit);
      is.open(file);
      int           c;
      // Display unchanged up to first header
      for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
        ++n_line;
        std::getline(is, line);
        if(args.line_number_flag)
          std::cout << n_line << ':';
        std::cout << line << '\n';
      }

      while(c != EOF && max_count > 0) {
        const size_t header_n_line = ++n_line;
        std::getline(is, header);
        bool header_printed = false;

        for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
          ++n_line;
          std::getline(is, line);
          if(regex_search(line.cbegin(), line.cend(), regexp) ^ args.invert_match_flag) {
            if(!header_printed) {
              if(args.line_number_flag)
                std::cout << header_n_line << ':';
              std::cout << header << '\n';
              header_printed = true;
            }
            if(args.line_number_flag)
              std::cout << n_line << ':';
            std::cout << line << '\n';
            if(--max_count == 0) break;
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
