#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>

#include <mummer/fasta.hpp>

// Return the reverse complement of sequence. This allows searching
// the plus strand of instances on the minus strand.
void reverse_complement(std::string &seq_rc, bool nucleotides_only) {
  // Reverse in-place.
  reverse(seq_rc.begin(), seq_rc.end());
  for(long i = 0; i < (long)seq_rc.length(); i++) {
    // Adapted from Kurtz code in MUMmer v3. 
    switch(seq_rc[i]) {
    case 'a': seq_rc[i] = 't'; break;
    case 'c': seq_rc[i] = 'g'; break;
    case 'g': seq_rc[i] = 'c'; break;
    case 't': seq_rc[i] = 'a'; break;
    case 'r': seq_rc[i] = 'y'; break; /* a or g */
    case 'y': seq_rc[i] = 'r'; break; /* c or t */
    case 's': seq_rc[i] = 's'; break; /* c or g */
    case 'w': seq_rc[i] = 'w'; break; /* a or t */
    case 'm': seq_rc[i] = 'k'; break; /* a or c */
    case 'k': seq_rc[i] = 'm'; break; /* g or t */
    case 'b': seq_rc[i] = 'v'; break; /* c, g or t */
    case 'd': seq_rc[i] = 'h'; break; /* a, g or t */
    case 'h': seq_rc[i] = 'd'; break; /* a, c or t */
    case 'v': seq_rc[i] = 'b'; break; /* a, c or g */
    default:  
      if(!nucleotides_only) seq_rc[i] = 'n'; 
      break; /* anything */
    }
  }
}


// Trim a string, giving start and end in trimmed version.
// NOTE: Assumes line.length() > 0!!!!
void trim(std::string &line, long &start, long &end) {
  // Trim leading spaces.
  for( ; start < (long)line.length() && line[start] == ' '; ++start) ;

  // Trim trailing spaces.
  for( ; end > 0 && line[end] == ' '; --end) ;
}


// Concatenate new sequences to set, keep track of lengths.
// NOTE: Concatenation using the '`' character to separate strings!
void load_fasta(const std::string& filename, std::string &S, std::vector<std::string> &descr, std::vector<long> &startpos) {
  std::string meta, line;

  // Everything starts at zero.
  std::ifstream data(filename.c_str());

  if(!data.is_open()) { std::cerr << "unable to open " << filename << std::endl; exit(1); }
  int c = data.peek();
  if(c != '>') { std::cerr << "first character must be a '>', got '" << (char)c << "'" << std::endl; exit(1); }

  while(c != EOF) {
    std::getline(data, line); // Load one line at a time.
    assert(line.back() != '\r');

    // Read metadata
    const size_t start = line.find_first_not_of(" ", 1);
    if(start != std::string::npos) {
      const size_t end = line.find_first_of(" ", start);
      descr.push_back(line.substr(start, std::min(end, line.size()) - start));
    }

    // Read sequence
    if(!S.empty())
      S += '`';
    startpos.push_back(S.length());
    for(c = data.peek(); c != EOF && c != '>'; c = data.peek()) {
      std::getline(data, line);
      const size_t start = line.find_first_not_of(" ");
      const size_t end = std::min(line.size(), line.find_last_not_of(" "));
      for(size_t i = start; i <= end; ++i)
        S += std::tolower(line[i]);
    }
  }
}


