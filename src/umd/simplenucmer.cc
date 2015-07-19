#include <iostream>
#include <fstream>
#include <umd/nucmer.hpp>

std::string read_sequence(const char* file) {
  std::string res;
  std::ifstream is(file);
  std::string line;

  std::getline(is, line); // Ignore header
  while(std::getline(is, line))
    res += line;
  return res;
}

int main(int argc, char *argv[]) {
  if(argc != 3) {
    std::cerr << "Usage: " << argv[0] << " ref qry" << std::endl;
    exit(1);
  }

  std::string ref = read_sequence(argv[1]);
  std::string qry = read_sequence(argv[2]);

  const auto alignments = mummer::nucmer::match(ref.c_str(), qry.c_str());
  mummer::postnuc::printDeltaAlignments(alignments,
                                        "ref", ref.size(), "qry", qry.size(),
                                        std::cout);
  return 0;
}
