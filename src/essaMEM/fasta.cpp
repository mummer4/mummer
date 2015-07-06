#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>

#include "fasta.hpp"

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
  for(long i = start; i < (int)line.length(); i++) { 
    if(line[i] != ' ') { start = i; break; } 
  }
  // Trim trailing spaces.
  for(long i = line.length() - 1; i >= 0; i--) { 
    if(line[i] != ' ') { end = i; break; } 
    if(i == 0) break;
  }
}


// Concatenate new sequences to set, keep track of lengths.
// NOTE: Concatenation using the '`' character to separate strings!
void load_fasta(std::string filename, std::string &S, std::vector<std::string> &descr, std::vector<long> &startpos) {
  std::string meta, line;
  long length = 0;

  // Everything starts at zero.
  startpos.push_back(0);

  std::ifstream data(filename.c_str());

  if(!data.is_open()) { std::cerr << "unable to open " << filename << std::endl; exit(1); } 

  while(!data.eof()) {
    getline(data, line); // Load one line at a time.
    if(line.length() == 0) continue;
    if(line[line.length()-1]=='\r'){
        line.erase(--line.end());
        std::cerr << "shouldn't happen" << std::endl;
        exit(1);
    }

    long start = 0, end = line.length() - 1;

    // Meta tag line and start of a new sequence.
    if(line[0] == '>') {
      // Save previous sequence and meta data.
      if(length > 0) {
	descr.push_back(meta);
	S += '`'; // ` character used to separate strings
	startpos.push_back(S.length());
	//lengths.push_back(length+1);
      }
      // Reset parser state.
      start = 1; meta = ""; length = 0;
    }
    trim(line, start, end);
    // Collect meta data.
    if(line[0] == '>') {
      for(long i = start; i <= end; i++) { if(line[i] == ' ') break; meta += line[i]; }
    }
    else { // Collect sequence data.
      length += end - start + 1;
      for(long i = start; i <= end; i++) { 
	S += std::tolower(line[i]);
      }
    }
  }
  if(length > 0) {
    descr.push_back(meta);
  }  
  std::cerr << "# S.length=" << S.length() << std::endl;
  for(long i = 0; i < (long)descr.size(); i++) {
    std::cerr << "# " << descr[i] << " " << startpos[i] << std::endl;
  }
}


