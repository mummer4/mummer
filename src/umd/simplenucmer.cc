#include <iostream>
#include <fstream>
#include <umd/nucmer.hpp>

std::string read_sequence(const char* file, std::string& header) {
  std::string res;
  std::ifstream is(file);
  std::string line;

  std::getline(is, line);
  header = line.substr(1);
  while(std::getline(is, line))
    res += line;
  return res;
}

int main(int argc, char *argv[]) {
  if(argc != 3) {
    std::cerr << "Usage: " << argv[0] << " ref qry" << std::endl;
    exit(1);
  }

  std::string ref_header, qry_header;
  std::string ref = read_sequence(argv[1], ref_header);
  std::string qry = read_sequence(argv[2], qry_header);

  const auto alignments = mummer::nucmer::align_sequences(ref.c_str(), qry.c_str());
  mummer::postnuc::printDeltaAlignments(alignments,
                                        ref_header, ref.size(), qry_header, qry.size(),
                                        std::cout);
  return 0;
}
