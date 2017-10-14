#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "check_cigar_cmdline.hpp"

typedef std::map<std::string, std::string> name_seq_map;

name_seq_map read_fasta(const char* path) {
  name_seq_map res;
  std::ifstream is(path);

  int c = is.peek();
  if(c != '>')
    check_cigar_cmdline::error() << "Invalid fasta file '" << path << '\'';
  std::string line;
  while(c != EOF) {
    std::getline(is, line);
    std::string& seq = res[line.substr(1, line.find_first_of(" \t", 1) - 1)];
    for(c = is.peek(); c != EOF && c != '>'; c = is.peek()) {
      std::getline(is, line);
      seq += line;
    }
  }

  return res;
}

char comp(char c) {
  switch(c) {
  case 'a': case 'A': return 't';
  case 'c': case 'C': return 'g';
  case 'g': case 'G': return 'c';
  case 't': case 'T': return 'a';
  default: return 'n';
  }
}

std::string rev_comp(const std::string& seq) {
  std::string res;
  for(auto it = seq.rbegin(); it != seq.rend(); ++it)
    res += comp(*it);
  return res;
}

std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> res;
  size_t pi = 0;
  size_t ni = line.find_first_of("\t");
  while(true) {
    res.push_back(line.substr(pi, ni - pi));
    if(ni == std::string::npos) break;
    pi = line.find_first_not_of("\t", ni);
    if(pi == std::string::npos) break;
    ni = line.find_first_of("\t", pi);
  }
  return res;
}

const std::string& find(const name_seq_map& map, const std::string name) {
    const auto it = map.find(name);
    if(it == map.end())
      check_cigar_cmdline::error() << "Did not sequence named '" << name << '\'';
    return it->second;
}

struct cig_info {
  const int  l;
  const char t;
};
std::vector<cig_info> split_cigar(const std::string& s) {
  std::vector<cig_info> res;
  size_t pi = 0;
  size_t ni = s.find_first_of("MIDNSHP=X");
  while(ni != std::string::npos) {
    const int l = std::stoi(s.substr(pi, ni - pi));
    res.push_back({ l, s[ni] });
    pi = ni + 1;
    ni = s.find_first_of("MIDNSHP=X", pi);
  }
  return res;
}

long find_nm(const std::vector<std::string>& fields) {
  for(unsigned i = 11; i < fields.size(); ++i) {
    if(fields[i].substr(0, 5) == "NM:i:") {
      return std::stol(fields[i].substr(5));
    }
  }
  return -1;
}

// Read a SAM file and two fasta files. It checks that the CIGAR
// string represent the difference between the sequences. Every pair
// of sequences mentioned in the SAM file must be present in the
// fasta files.
int main(int argc, char *argv[]) {
  check_cigar_cmdline args(argc, argv);
  auto ref_map = read_fasta(args.ref_arg);
  auto qry_map = read_fasta(args.qry_arg);

  std::ifstream is(args.sam_arg);
  if(!is.good())
    check_cigar_cmdline::error() << "Invalid SAM file '" << args.sam_arg << '\'';
  std::string line;

  for(size_t linenb = 1; std::getline(is, line); ++linenb) {
    if(line[0] == '@') continue;
    const auto fields = split(line);

    const auto&       ref   = find(ref_map, fields[2]);
    const auto&       qry   = find(qry_map, fields[0]);
    const bool        rev   = std::stoul(fields[1]) & 0x10;
    const std::string sam_seq = fields[9];
    const std::string seq   = rev ? rev_comp(qry) : qry;
    auto              cigar = split_cigar(fields[5]);
    const long        pos   = std::stol(fields[3]);
    const long        nm    = find_nm(fields);

    if(nm < 0)
      check_cigar_cmdline::error() << "linenb:" << linenb << " did not find NM field";

    if(pos > ref.size() || pos < 1)
      check_cigar_cmdline::error() << "linenb:" << linenb << " pos " << pos << " invalid for ref " << fields[2]
                                   << " of size " << ref.size();

    long refp  = pos - 1;
    long qryp  = 0;
    long sqryp = 0;
    long diffs = 0;
    for(const auto& ci : cigar) {
      if(refp > ref.size())
        check_cigar_cmdline::error() << "linenb:" << linenb << " pos " << refp << " beyond end of reference " << fields[2]
                                     << " of size " << ref.size();
      if(qryp > seq.size())
        check_cigar_cmdline::error() << "linenb:" << linenb << " pos " << qryp << " beyond end of qry " << fields[0]
                                     << " of size " << ref.size();

      switch(ci.t) {
      case 'H':
        qryp += ci.l;
        continue;
      case 'S':
        qryp  += ci.l;
        sqryp += ci.l;
        continue;
      case 'I':
        qryp  += ci.l;
        sqryp += ci.l;
        diffs += ci.l;
        continue;
      case 'D':
        refp  += ci.l;
        diffs += ci.l;
        continue;
      case 'M':
        break;
      default:
        check_cigar_cmdline::error() << "linenb:" << linenb << " unsupported type " << ci.t;
        break;
      }
      if(refp + ci.l > ref.size())
        check_cigar_cmdline::error() << "linenb:" << linenb << " match at " << refp << ':' << ci.l
                                     << " goes beyond end of reference " << fields[2]
                                     << " of size " << ref.size();
      if(qryp + ci.l > seq.size())
        check_cigar_cmdline::error() << "linenb:" << linenb << " match at " << qryp << ':' << ci.l
                                     << " goes beyond end of qry " << fields[0]
                                     << " of size " << seq.size();
      if(sam_seq[0] != '*' && sqryp >= sam_seq.size())
        check_cigar_cmdline::error() << "linenb:" << linenb << " match at " << qryp << ':' << ci.l
                                     << " goes beyond end of sequence in sam file " << sam_seq
                                     << " of size " << sam_seq.size();
      for(int i = 0; i < ci.l; ++i, ++refp, ++qryp, ++sqryp) {
        diffs += std::tolower(ref[refp]) != std::tolower(seq[qryp]);
        if(sam_seq[0] != '*' && std::tolower(qry[qryp]) != std::tolower(sam_seq[sqryp]))
            check_cigar_cmdline::error() << "linenb:" << linenb << " qry sequence does not match sequence in sam file at pos "
                                         << qryp << " and " << sqryp;
      }
    }
    if(diffs != nm)
      check_cigar_cmdline::error() << "linenb:" << linenb << " invalid NM field " << nm << ", should be " << diffs;
  }

  return 0;
}
