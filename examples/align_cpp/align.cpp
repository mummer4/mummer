#include <iostream>
#include <fstream>
#include <mummer/nucmer.hpp>

std::string read_file(const char* path) {
  std::ifstream is(path);
  std::string res, line;
  while(std::getline(is, line))
    res += line;
  return res;
}

int main(int argc, char *argv[]) {
  std::string ref = read_file(argv[1]);
  std::string qry = read_file(argv[2]);

  mummer::nucmer::Options o;
  o.minmatch(10).mincluster(10);
  auto aligns = mummer::nucmer::align_sequences(ref, qry, o);
  for(const auto& a : aligns) {
    std::cout << a.sA << ' ' << a.eA << ' ' << a.sB << ' ' << a.eB << ' '
              << a.Errors << ' ' << a.SimErrors << ' ' << a.NonAlphas << '\n';
    for(auto d : a.delta)
      std::cout << d << '\n';
    std::cout << "0\n";
  }

  return 0;
}
