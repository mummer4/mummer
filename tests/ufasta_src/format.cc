#include <cstdlib>
#include <string>
#include <fstream>
#include <cctype>

#include <format_cmdline.hpp>

void cleanup(std::string& l, const bool lc, const bool uc, const bool space) {
  if(!(lc || uc || space)) return;
  size_t j = 0;
  for(size_t i = 0; i < l.size(); ++i) {
    if(space && std::isspace(l[i])) continue;
    if(lc) l[j] = std::tolower(l[i]);
    if(uc) l[j] = std::toupper(l[i]);
    ++j;
  }
  l.resize(j);
}

int format_main(int argc, char *argv[]) {
  format_cmdline args(argc, argv);
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

      int c;
      // Display unchanged up to first header
      for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
        std::getline(is, line);
        std::cout << line << '\n';
      }

      while(c != EOF) {
        std::getline(is, line); // Header
        std::cout << line << '\n';

        size_t ooff  = 0;       // current out offset
        bool   empty = true;
        for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
          std::getline(is, line);
          cleanup(line, args.lower_case_flag, args.upper_case_flag, args.spaces_flag);
          size_t loff = 0; // current line offset
          while(true) {
            const size_t out_left  = args.line_length_arg - ooff;
            const size_t line_left = line.size() - loff;
            if(out_left < line_left) {
              std::cout.write(line.c_str() + loff, out_left);
              std::cout << '\n';
              empty = false;
              ooff  = 0;
              loff += out_left;
            } else {
              std::cout.write(line.c_str() + loff, line_left);
              empty = false;
              ooff += line_left;
              break;
            }
          }
        }
        if(!empty)
          std::cout << '\n';
      }
    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
