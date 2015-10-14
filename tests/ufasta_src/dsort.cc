#include <fstream>
#include <algorithm>

#include <common.hpp>
#include <dsort_cmdline.hpp>

int dsort_main(int argc, char *argv[]) {
  dsort_cmdline args(argc, argv);
  if(args.file_arg.empty()) {
    args.file_arg.push_back("/dev/stdin");
    if(isatty(0))
      std::cerr << "Warning: reading from terminal" << std::endl;
  }

  int         res = EXIT_SUCCESS;
  int         c;
  entry       data;
  std::string header;
  for(const auto& file : args.file_arg) {
    try {
      std::ifstream is;
      is.exceptions(std::ios::failbit|std::ios::badbit);
      is.open(file);
      // Display unchanged up to first header
      for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
        std::getline(is, header);
        std::cout << header << '\n';
      }

      while(c != EOF) {
        std::getline(is, header);
        data.size = 0;
        for(c = is.peek() ; c != '>' && c != EOF; c = is.peek())
          data.add_line(is);
        std::sort(data.lines.begin(), data.lines.begin() + data.size);
        std::cout << header << '\n';
        for(int i = 0; i < data.size; ++i)
          std::cout << data.lines[i] << '\n';
      }
    } catch(std::ios::failure) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      res = EXIT_FAILURE;
    }
  }

  return res;
}
