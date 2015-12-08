#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <signal.h>
#include <unistd.h>

#include <cstdlib>
#include <cstring>
#include <iostream>

typedef int (main_func_t)(int argc, char *argv[]);

main_func_t one_main;
// main_func_t sizes_main;
// main_func_t head_main;
// main_func_t tail_main;
// main_func_t rc_main;
// main_func_t n50_main;
// main_func_t extract_main;
// main_func_t format_main;
main_func_t sort_main;
main_func_t dsort_main;

// #ifdef HAVE_BOOST_REGEX
// main_func_t hgrep_main;
// main_func_t dgrep_main;
// #endif

main_func_t sos;
main_func_t version;

struct cmd_func {
  const char  *cmd;
  main_func_t *func;
};
cmd_func cmd_list[] = {
  // {"one",               &one_main},
  // {"sizes",             &sizes_main},
  // {"head",              &head_main},
  // {"tail",              &tail_main},
  // {"rc",                &rc_main},
  // {"n50",               &n50_main},
  // {"extract",           &extract_main},
  // {"format",            &format_main},
  {"sort",              &sort_main},
  {"hsort",             &sort_main},
  {"dsort",             &dsort_main},
// #ifdef HAVE_BOOST_REGEX
//   {"hgrep",             &hgrep_main},
//   {"dgrep",             &dgrep_main},
// #endif

  /* help in all its form. Must be first non-command */
  {"help",              &sos},
  {"-h",                &sos},
  {"-help",             &sos},
  {"--help",            &sos},
  {"-?",                &sos},
  {"--version",         &version},
  {"-V",                &version},
  {"",                  0}
};



void __sos(std::ostream *os)
{
  *os << "Usage: ufasta <cmd> [options] arg..."  << std::endl <<
    "Where <cmd> is one of: ";
  bool comma = false;
  for(cmd_func *ccmd = cmd_list; ccmd->func != sos; ccmd++) {
    *os << (comma ? ", " : "") << ccmd->cmd;
    comma = true;
  }
  *os << "." << std::endl;
  *os << "Options:" << std::endl <<
    "  --version        Display version" << std::endl <<
    "  --help           Display this message" << std::endl;
}

int sos(int argc, char *argv[])
{
  __sos(&std::cout);
  return 0;
}

int version(int argc, char *argv[])
{
#ifdef PACKAGE_STRING
  std::cout << PACKAGE_STRING << std::endl;
#else
  std::cout << "0.0.0" << std::endl;
#endif
  return 0;
}

void sigpipe_handler(int sig) {
  _exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
  std::string error;

  // Ignore SIGPIPE. It causes ufasta to fail if output is sent to a
  // pipe, which is not very useful for us. Simply exit successfully.
  {
    struct sigaction sig;
    memset(&sig, '\0', sizeof(sig));
    sig.sa_handler = sigpipe_handler;
    if(sigaction(SIGPIPE, &sig, nullptr) == -1)
      perror("sigaction");
  }

  if(argc < 2) {
    error = "Too few arguments";
  } else {
    for(cmd_func *ccmd = cmd_list; ccmd->func != 0; ccmd++) {
      if(!strcmp(ccmd->cmd, argv[1]))
        return ccmd->func(argc - 1, argv + 1);
    }
    error = "Unknown command '";
    error += argv[1];
    error += "'\n";
  }

  std::cerr << error << std::endl;
  __sos(&std::cerr);
  return EXIT_FAILURE;
}
