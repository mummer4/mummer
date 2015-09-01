#include <iostream>
#include <fstream>
#include <vector>
#include <mummer/nucmer.hpp>

std::string read_sequence(std::istream& is, std::string& header) {
  std::string res;
  std::string line;

  std::getline(is, line);
  header = line.substr(1);
  while(is.peek() != '>' && std::getline(is, line))
    res += line;
  return res;
}

int main(int argc, char *argv[]) {
  if(argc != 3) {
    std::cerr << "Usage: " << argv[0] << " ref qry" << std::endl;
    exit(1);
  }

  std::vector<std::string> qry_sequences;
  std::vector<std::string> qry_headers;
  { std::ifstream qrys(argv[2]);
    while(qrys.peek() == '>') {
      std::string header;
      std::string seq = read_sequence(qrys, header);
      qry_sequences.push_back(std::move(seq));
      qry_headers.push_back(std::move(header));
    }
  }

  std::ifstream refs(argv[1]);
  while(refs.peek() == '>') {
    std::string header;
    const std::string seq = read_sequence(refs, header);
    // std::cout << '>' << header << ' ' << seq.size() << '\n'
    //           << seq << '\n';
    mummer::nucmer::SequenceAligner aligner(seq.c_str());
    for(size_t i = 0; i < qry_sequences.size(); ++i) {
      std::cout << i << std::endl;
      mummer::postnuc::printDeltaAlignments(aligner.align(qry_sequences[i].c_str(), qry_sequences[i].size()),
                                            header, seq.size(), qry_headers[i], qry_sequences[i].size(),
                                            std::cout);
      std::cout << std::flush;
    }
  }

  return 0;
}
