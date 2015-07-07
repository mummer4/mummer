#ifndef __FASTA_HPP__
#define __FASTA_HPP__

#include <string>
#include <vector>

// using namespace std;

void reverse_complement(std::string &seq_rc, bool nucleotides_only);
void trim(std::string &line, long &start, long &end);
void load_fasta(const std::string& filename, std::string &S, std::vector<std::string> &descr, std::vector<long> &startpos);

#endif // __FASTA_HPP__

