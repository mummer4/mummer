#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <cerrno>
#include <rsort_cmdline.hpp>
#include "common.hpp"

struct sequence {
  ssize_t index;
  std::vector<std::string> data;
  void clear() { data.clear(); }
};
std::ostream& operator<<(std::ostream& os, const sequence& seq) {
  for(const auto& s : seq.data)
    os << s << '\n';
  return os;
}
// Less operator for the heap. It sorts indices in reverse order
// because we want a min-heap, not a max-heap
bool operator<(const sequence& s1, const sequence& s2) {
  return s1.index > s2.index;
}

int rsort_main(int argc, char *argv[]) {
  rsort_cmdline args(argc, argv);

  // Heap of sequence information
  size_t heap_size = 0;
  ssize_t cindex = args.start_arg; // Index of next sequence to output
  std::vector<sequence> sequences;

  for(const auto path : args.file_arg) {
    std::ifstream is(path);
    if(!is.good())
      rsort_cmdline::error() << "Failed to open file '" << path << '\'';

    // Skip to first header
    int c;
    for(c = is.peek(); c != '>' && c != EOF; c = is.peek())
      skip_line(is);

    while(c != EOF) {
      if(heap_size == sequences.size())
        sequences.resize(heap_size + 1);
      sequence& cseq = sequences[heap_size];
      cseq.clear();

      cseq.data.resize(1);
      std::getline(is, cseq.data[0]);
      const char* key = find_token(args.key_arg - 1, cseq.data.front().c_str());
      if(!key)
        rsort_cmdline::error() << "Sequence header without a key: " << cseq.data.front();
      errno = 0;
      if(args.prefix_flag) {
        while(*key && !isspace(*key) && !isdigit(*key)) ++key;
        if(!isdigit(*key))
          rsort_cmdline::error() << "No number after key prefix:: " << cseq.data.front();
      }
      ssize_t index = strtoll(key, nullptr, 10);
      if(errno)
        rsort_cmdline::error() << "Key is not a number: " << cseq.data.front();
      if(index < cindex)
        rsort_cmdline::error() << "Repeated key: " << cseq.data.front();
      cseq.index = index;

      for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
        cseq.data.resize(cseq.data.size() + 1);
        std::getline(is, cseq.data.back());
      }

      // Decide what to do with the sequence. If the next one to
      // output, do so and output everything that should come next
      // from the heap. Otherwise, send through the heap.
      if(index == cindex) {
        std::cout << cseq;
        for(++cindex; heap_size && sequences.front().index == cindex; ++cindex, --heap_size) {
          std::cout << sequences.front();
          std::pop_heap(sequences.begin(), sequences.begin() + heap_size);
        }
      } else {
        ++heap_size;
        std::push_heap(sequences.begin(), sequences.begin() + heap_size);
      }
    }
  }
  if(heap_size) // Heap not empty: we are missing indices
    rsort_cmdline::error() << "Missing key: " << cindex;

  return EXIT_SUCCESS;
}
