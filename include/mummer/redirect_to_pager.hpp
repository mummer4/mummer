#ifndef __REDIRECT_TO_PAGER_H__
#define __REDIRECT_TO_PAGER_H__

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Class that will redirect stdout to a pager if: (1) stdout is a tty
// and (2) it can launch a pager. The pager is the content of the
// PAGER environment variables or "more" if not set. This works with
// stdio.
struct stdio_launch_pager {
  static constexpr const char* dflt_pager = "more";
  FILE* pager_handle;

  stdio_launch_pager()
    : pager_handle(NULL)
  {
    if(!isatty(1)) return;
    // Open a pager
    const char* pager      = getenv("PAGER");
    if(pager)
      pager_handle = start_pager(pager);
    if(!pager_handle) {
      pager_handle = start_pager(dflt_pager);
      if(!pager_handle)
        return;
    }

    // Close stdout and put pager in its place
    if(fileno(pager_handle) == 1)
      return; // Done already!
    if(dup2(fileno(pager_handle), 1) == -1)
      stop_pager();
  }

  ~stdio_launch_pager() { stop_pager(); }

  FILE* start_pager(const char* cmd) {
    return popen(cmd, "w");
  }

  void stop_pager() {
    if(pager_handle) {
      pclose(pager_handle);
      pager_handle = NULL;
    }
  }
};

#endif /* __REDIRECT_TO_PAGER_H__ */
