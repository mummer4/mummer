#include "common.hpp"
#include <split_cmdline.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/select.h>


struct output_pipe {
  int               fd;
  size_t            pos, ipos; // pos to write from buffer and ipos to read in buffer
  std::vector<char> buffer;

  output_pipe() : fd(-1), pos(0) { }
  output_pipe(const char* path)
    : fd(open(path, O_WRONLY|O_CREAT|O_NONBLOCK|O_APPEND, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))
    , pos(0)
    , ipos(0)
    , buffer(1024, '\0')
  { }
  output_pipe(output_pipe&& rhs)
    : fd(rhs.fd)
    , pos(rhs.pos)
    , ipos(rhs.ipos)
    , buffer(std::move(rhs.buffer))
  {
    rhs.fd = -1;
  }
  ~output_pipe() {
    close();
  }

  output_pipe& operator=(output_pipe&& rhs) {
    fd     = rhs.fd;
    rhs.fd = -1;
    pos    = rhs.pos;
    ipos   = rhs.ipos;
    buffer = std::move(rhs.buffer);
    return *this;
  }

  void close() {
    if(fd >= 0)
      ::close(fd);
    fd = -1;
  }

  bool append_sequence(std::istream& is) {
    if(pos < ipos) { // Try to write if have stuff remaining in buffer
      while(true) {
        ssize_t res = write(fd, buffer.data() + pos, ipos - pos);
        if(res == -1) {
          if(errno == EAGAIN || errno == EWOULDBLOCK) return true;
          if(errno == EINTR) continue;
          std::cerr << "Warning: error writing to pipe: " << strerror(errno);
          return false;
        }
        pos += res;
        break;
      }
      if(pos < ipos) return true;
    }

    // Need to refill buffer
    pos = 0;
    ipos = append_line(is, buffer, 0);
    while(is.good() && is.peek() != EOF && is.peek() != '>')
      ipos = append_line(is, buffer, ipos);
    return ipos > 0;
  }
};

int split_main(int argc, char* argv[]) {
  split_cmdline args(argc, argv);

  std::vector<output_pipe> pipes;
  if(args.output_arg.empty())
    args.output_arg.push_back("/dev/stdout");
  for(auto& path : args.output_arg) {
    output_pipe pipe(path);
    if(pipe.fd == -1)
      split_cmdline::error() << "Failed to open file '" << path << "':" << strerror(errno);
    pipes.push_back(std::move(pipe));
  }

  if(!strcmp(args.input_arg, "/dev/stdin") && isatty(0))
      std::cerr << "Warning: reading from terminal" << std::endl;

  const char* input_path = args.input_arg;
  std::ifstream input(input_path);
  if(!input.good())
    split_cmdline::error() << "Failed to open input file '" << input_path << "'";
  std::string buffer;

  while(true) {
    fd_set output_set;
    FD_ZERO(&output_set);
    int fd_max = -1;
    for(auto& pipe : pipes) {
      if(pipe.fd != -1)
        FD_SET(pipe.fd, &output_set);
      fd_max = std::max(fd_max, pipe.fd);
    }
    if(fd_max == -1)
      break;
    ++fd_max;

    while(true) {
      int res = select(fd_max, nullptr, &output_set, nullptr, nullptr);
      if(res == -1 && errno == EINTR) continue;
      if(res == -1)
        split_cmdline::error() << "Error while doing select: " << strerror(errno);
      break;
    }

    for(auto& pipe : pipes) {
      if(pipe.fd == -1) continue;
      if(FD_ISSET(pipe.fd, &output_set) && !pipe.append_sequence(input))
        pipe.close();
    }
  }

  return EXIT_SUCCESS;
}
