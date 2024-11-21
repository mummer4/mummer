#include <cstdlib>
#include <string>
#include <fstream>

#include <one_cmdline.hpp>

int one_main(int argc, char *argv[]) {
  one_cmdline args(argc, argv);
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

        for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
          std::getline(is, line);
          std::cout << line;
        }
        std::cout << '\n';
      }
    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
