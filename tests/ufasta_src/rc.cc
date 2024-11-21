#include <cstdlib>
#include <string>
#include <fstream>
#include <cctype>

#include <rc_cmdline.hpp>

const char rc[256] = {
 (char)  0, (char)  1, (char)  2, (char)  3, (char)  4, (char)  5, (char)  6, (char)  7,
 (char)  8, (char)  9, (char) 10, (char) 11, (char) 12, (char) 13, (char) 14, (char) 15,
 (char) 16, (char) 17, (char) 18, (char) 19, (char) 20, (char) 21, (char) 22, (char) 23,
 (char) 24, (char) 25, (char) 26, (char) 27, (char) 28, (char) 29, (char) 30, (char) 31,
 (char) 32, (char) 33, (char) 34, (char) 35, (char) 36, (char) 37, (char) 38, (char) 39,
 (char) 40, (char) 41, (char) 42, (char) 43, (char) 44, (char) 45, (char) 46, (char) 47,
 (char) 48, (char) 49, (char) 50, (char) 51, (char) 52, (char) 53, (char) 54, (char) 55,
 (char) 56, (char) 57, (char) 58, (char) 59, (char) 60, (char) 61, (char) 62, (char) 63,
 (char) 64, (char) 84, (char) 66, (char) 71, (char) 68, (char) 69, (char) 70, (char) 67,
 (char) 72, (char) 73, (char) 74, (char) 75, (char) 76, (char) 77, (char) 78, (char) 79,
 (char) 80, (char) 81, (char) 82, (char) 83, (char) 65, (char) 85, (char) 86, (char) 87,
 (char) 88, (char) 89, (char) 90, (char) 91, (char) 92, (char) 93, (char) 94, (char) 95,
 (char) 96, (char)116, (char) 98, (char)103, (char)100, (char)101, (char)102, (char) 99,
 (char)104, (char)105, (char)106, (char)107, (char)108, (char)109, (char)110, (char)111,
 (char)112, (char)113, (char)114, (char)115, (char) 97, (char)117, (char)118, (char)119,
 (char)120, (char)121, (char)122, (char)123, (char)124, (char)125, (char)126, (char)127,
 (char)128, (char)129, (char)130, (char)131, (char)132, (char)133, (char)134, (char)135,
 (char)136, (char)137, (char)138, (char)139, (char)140, (char)141, (char)142, (char)143,
 (char)144, (char)145, (char)146, (char)147, (char)148, (char)149, (char)150, (char)151,
 (char)152, (char)153, (char)154, (char)155, (char)156, (char)157, (char)158, (char)159,
 (char)160, (char)161, (char)162, (char)163, (char)164, (char)165, (char)166, (char)167,
 (char)168, (char)169, (char)170, (char)171, (char)172, (char)173, (char)174, (char)175,
 (char)176, (char)177, (char)178, (char)179, (char)180, (char)181, (char)182, (char)183,
 (char)184, (char)185, (char)186, (char)187, (char)188, (char)189, (char)190, (char)191,
 (char)192, (char)193, (char)194, (char)195, (char)196, (char)197, (char)198, (char)199,
 (char)200, (char)201, (char)202, (char)203, (char)204, (char)205, (char)206, (char)207,
 (char)208, (char)209, (char)210, (char)211, (char)212, (char)213, (char)214, (char)215,
 (char)216, (char)217, (char)218, (char)219, (char)220, (char)221, (char)222, (char)223,
 (char)224, (char)225, (char)226, (char)227, (char)228, (char)229, (char)230, (char)231,
 (char)232, (char)233, (char)234, (char)235, (char)236, (char)237, (char)238, (char)239,
 (char)240, (char)241, (char)242, (char)243, (char)244, (char)245, (char)246, (char)247,
 (char)248, (char)249, (char)250, (char)251, (char)252, (char)253, (char)254, (char)255
};

void rc_sequence(std::string& seq) {
  ssize_t i = 0, j = (ssize_t)seq.size() - 1;
  for( ; i < j; ++i, --j) {
    char tmp = seq[i];
    seq[i] = rc[(int)seq[j]];
    seq[j] = rc[(int)tmp];
  }
  if(i == j)
    seq[i] = rc[(int)seq[i]];
}

bool is_canonical(const std::vector<std::string>& sequences, size_t nb_lines) {
  if(nb_lines == 0) return true; // vacuously canonical

  auto vst = sequences.cbegin();
  auto ven = vst;
  std::advance(ven, nb_lines - 1);

  auto st = vst->cbegin();
  auto en = ven->cend();

  while(true) {
    while(st == vst->cend()) {
      ++vst; // guaranteed vst < sequences.cend()
      st = vst->cbegin();
    }
    while(en == ven->cbegin()) {
      --ven; // guaranteed ven >= sequences.cbegin()
      en = ven->cend();
    }

    --en;
    if(vst >= ven || (vst == ven && st >= en)) break;
    const char c1 = std::tolower(*st);
    const char c2 = std::tolower(rc[(int)*en]);
    ++st;

    if(c1 < c2) return true;
    if(c1 > c2) return false;
  }

  return true; // palyndromic!
}

int rc_main(int argc, char *argv[]) {
  rc_cmdline args(argc, argv);
  if(args.file_arg.empty()) {
    args.file_arg.push_back("/dev/stdin");
    if(isatty(0))
      std::cerr << "Warning: reading from terminal" << std::endl;
  }


  std::string header;
  std::vector<std::string> sequences;
  for(auto file : args.file_arg) {
    try {
      std::ifstream is;
      is.exceptions(std::ios::failbit|std::ios::badbit);
      is.open(file);

      int c;
      // Display unchanged up to first header
      for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
        std::getline(is, header);
        std::cout << header << '\n';
      }

      while(c != EOF) {
        std::getline(is, header);
        std::cout << header << '\n';

        size_t nb_lines = 0;
        for(c = is.peek(); c != '>' && c != EOF; c = is.peek(), ++nb_lines) {
          if(nb_lines >= sequences.size())
            sequences.push_back("");
          std::getline(is, sequences[nb_lines]);
        }

        if(!args.canonical_flag || !is_canonical(sequences, nb_lines)) {
          for(size_t i = nb_lines; i > 0; --i) {
            rc_sequence(sequences[i - 1]);
            std::cout << sequences[i - 1] << '\n';
          }
        } else {
          for(size_t i = 0; i < nb_lines; ++i)
            std::cout << sequences[i] << '\n';
        }
      }
    } catch(std::ios::failure&) {
      std::cerr << "Error with file '" << file << '\'' << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
