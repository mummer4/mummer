#include <iostream>
#include <fstream>
#include <mummer/nucmer.hpp>

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

  std::string qry_header;
  std::string qry = read_sequence(argv[2], qry_header);

  mummer::nucmer::FileAligner aligner(argv[1], mummer::nucmer::Options().mincluster(20).minmatch(10));
  aligner.align(qry.c_str(),
                [&](std::vector<mummer::postnuc::Alignment>&& als,
                   const mummer::nucmer::FastaRecordPtr& Af, const mummer::nucmer::FastaRecordSeq & Bf) {
                  //  const auto alignments = mummer::nucmer::align_sequences(ref.c_str(), qry.c_str());
                  mummer::postnuc::printDeltaAlignments(als, Af.Id(), Af.len(), qry_header, qry.size(),
                                                        std::cout);
                });
  return 0;
}
