#include <string>
#include <fstream>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <random>

#include "common.hpp"
#include <extract_cmdline.hpp>


std::unordered_set<std::string> read_names(const extract_cmdline& args) {
  std::unordered_set<std::string> res;

  for(const auto& n : args.name_arg)
    res.insert(n);
  std::string entry;
  for(const auto& file : args.name_file_arg) {
    std::ifstream is(file);
    if(is.fail())
      extract_cmdline::error() << "Error opening file '" << file << "' " << is.bad();
    for(is >> entry; is.good(); is >> entry)
      res.insert(entry);
  }
  return res;
}

double check_proba(const extract_cmdline& args) {
  if(!args.probabilistic_given) return 0.0;
  if(args.probabilistic_arg <= 0.0 || args.probabilistic_arg >= 1.0)
    extract_cmdline::error() << "Invalid probability " << args.probabilistic_arg << ", it must be in (0.0, 1.0)";
  return args.probabilistic_arg;
}

int extract_main(int argc, char *argv[]) {
  extract_cmdline args(argc, argv);
  if(args.file_arg.empty()) {
    args.file_arg.push_back("/dev/stdin");
    if(isatty(0))
      std::cerr << "Warning: reading from terminal" << std::endl;
  }

  const auto entries    = read_names(args);
  const auto end        = entries.cend();
  const bool no_entries = entries.empty();

  size_t max_count = args.max_count_given ? args.max_count_arg : std::numeric_limits<size_t>::max();
  if(max_count == 0) return EXIT_SUCCESS;

  std::default_random_engine             generator(time(nullptr)); // Weak random iniatilization. Ok for our purpose
  std::uniform_real_distribution<double> dist; // uniform in [0.0, 1.0)
  const double                           proba = check_proba(args);

  const bool range_limit = args.start_given || args.end_given;
  const size_t range_start = args.start_given ? args.start_arg : 0;
  const size_t range_end = args.end_given ? args.end_arg : std::numeric_limits<size_t>::max();
  if(range_start >= range_end)
    extract_cmdline::error() << "Range [" << range_start << ", " << range_end << ") specification error. Start of range must be less than end of range";

  // std::cerr << no_entries << ' ' << (!proba) << ' '
  //           << args.invert_match_flag << ' ' << range_limit << std::endl;

  if(no_entries && !proba) { // Nothing specific to extract
    if(!args.invert_match_flag) // Nothing to extract at all
      return EXIT_SUCCESS;

    if(!range_limit) {
      // Extract everything -> devolve to /bin/cat
      args.file_arg.insert(args.file_arg.begin(), "cat");
      args.file_arg.push_back(nullptr);
      execvp("cat", (char**)args.file_arg.data());
      return EXIT_FAILURE;
    }
  }

  std::string line;
  for(auto file : args.file_arg) {
    try {
      std::ifstream is;
      is.exceptions(std::ios::failbit|std::ios::badbit);
      is.open(file);
      int           c;
      // Display unchanged up to first header
      for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
        std::getline(is, line);
        std::cout << line << '\n';
      }

      while(c != EOF) {
        std::getline(is, line);
        const auto space = line.find_first_of(" \t", 1);
        const bool skip = ((!no_entries && entries.find(line.substr(1, space - 1)) == end) ||
                           (proba && dist(generator) >= proba)) ^ args.invert_match_flag;
        if(skip) {
          for(c = is.peek(); c != '>' && c != EOF; c = is.peek())
            skip_line(is);
          continue;
        }

        std::cout << line << '\n';
        if(!range_limit) {
          for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
            std::getline(is, line);
            std::cout << line << '\n';
          }
        } else {
          size_t byte = 0;
          for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
            std::getline(is, line);
            const size_t nbyte = byte + line.size();
            if(nbyte >= range_start && byte < range_end) { // need to display something
              const size_t s = range_start <= byte ? 0 : range_start - byte;
              const size_t l = nbyte < range_end ? std::string::npos : (range_end - std::max(range_start, byte));
              if(l)
                std::cout << line.substr(s, l) << '\n';
            }
            byte = nbyte;
          }
        }

        if(--max_count == 0) break;
      }
    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
