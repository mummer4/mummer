#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#include <type_traits>
#include <string>
#include <algorithm>

#include <sort_cmdline.hpp>

static sort_cmdline args;

struct close_fd {
  int         fd;
  char*       ptr;
  off_t       size;
  enum mem_type { MALLOC, MMAP };
  mem_type    type;

  close_fd() : fd(-1), ptr(nullptr), size(0) { }
  close_fd(int i) : fd(i), ptr(nullptr), size(0) { }
  close_fd(const close_fd&) = delete;
  close_fd(close_fd&& rhs) noexcept : fd(rhs.fd), ptr(rhs.ptr), size(rhs.size) {
    rhs.fd   = -1;
    rhs.ptr  = nullptr;
    rhs.size = 0;
  }
  ~close_fd() {
    close();
    unmap();
  }
  void close() {
    if(fd >= 0) ::close(fd);
    fd = -1;
  }
  void unmap() {
    if(ptr != nullptr) munmap((void*)ptr, size);
  }
};

template<typename T>
struct header_type {
  T           value;
  const char* start;
  const char* end;
  header_type(T v, const char* s, const char* e) : value(v), start(s), end(e) { }
};

#ifdef MREMAP_MAYMOVE
void* mem_alloc(size_t size) {
  return mmap(0, size, PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
}
void* mem_realloc(void* ptr, size_t old_size, size_t new_size) {
  return mremap(ptr, old_size, new_size, MREMAP_MAYMOVE);
}
#else
void* mem_alloc(size_t size) {
  return malloc(size);
}
void* mem_realloc(void* ptr, size_t old_size, size_t new_size) {
  return realloc(ptr, new_size);
}
#endif

void slurp_in(close_fd& fd) {
  fd.size = 1024 * 1024;
  fd.ptr  = (char*)mem_alloc(fd.size);
  if(fd.ptr == MAP_FAILED) return;

  size_t offset = 0;
  while(true) {
    size_t left = fd.size - offset - 1;
    while(left > 0) {
      ssize_t res = read(fd.fd, fd.ptr + offset, left);
      if(res == -1 && errno == EINTR) continue;
      if(res == -1) goto error;
      if(res == 0) {
        size_t new_size = offset + 1;
        void*  new_ptr  = mem_realloc(fd.ptr, fd.size, new_size);
        if(new_ptr == MAP_FAILED) goto error;
        fd.size = new_size;
        fd.ptr = (char*)new_ptr;
        return;
      }
      offset += res;
      left   -= res;
    }
    size_t new_size = fd.size * 2;
    void*  new_ptr  = mem_realloc(fd.ptr, fd.size, new_size);
    if(new_ptr == MAP_FAILED) goto error;
    fd.size = new_size;
    fd.ptr  = (char*)new_ptr;
  }

 error:
  int save_errno = errno;
  fd.unmap();
  errno = save_errno;
}

close_fd open_mmap(const char* path) {
  close_fd res(open(path, O_RDONLY));
  if(res.fd == -1)
    sort_cmdline::error() << "Failed to open file '" << path << "':" << strerror(errno);
  struct stat file_stat;
  memset(&file_stat, 0, sizeof(file_stat));
  if(fstat(res.fd, &file_stat) == -1)
    sort_cmdline::error() << "Failed to stat file '" << path << "':" << strerror(errno);
  res.size = file_stat.st_size + 1;
  res.ptr = (char*)mmap(nullptr, res.size, PROT_READ, MAP_PRIVATE, res.fd, 0);
  if(res.ptr == MAP_FAILED) {
    slurp_in(res);
    if(res.ptr == MAP_FAILED)
      sort_cmdline::error() << "Failed to read file '" << path << "':" << strerror(errno);
  }
  res.close();
  return res;
}

template<typename T>
struct header_traits { };

// Key is an signed int. Sort numerically
template<>
struct header_traits<int64_t> {
  typedef header_type<int64_t> type;

  static type name(const char* start, const char* end, const char* col) {
    static std::string buff;
    const char* const  ptr = std::min(col, end);
    const size_t       s   = strcspn(ptr, " \n\t");
    buff.assign(ptr, s);
    return type(std::stoll(buff), start, end);
  }

  static void sort(std::vector<type>& headers) {
    std::sort(headers.begin(), headers.end(), [](const type& x, const type& y) { return x.value < y.value; });
  }
};

// Key is a C string (type size_t is its length). Sort alphabetically.
struct str_type {
  const char* col;
  size_t      size;
};
template<>
struct header_traits<str_type> {
  typedef header_type<str_type> type;

  static type name(const char* start, const char* end, const char* col) {
    str_type s = { col, strcspn(std::min(col, end), args.header_full_flag ? "\n" : " \n\t") };
    return type(s, start, end);
  }

  static void sort(std::vector<type>& headers) {
    std::sort(headers.begin(), headers.end(),
              [](const type& x, const type& y) -> bool {
                const bool xshort = x.value.size < y.value.size;
                const int  res    = memcmp(x.value.col, y.value.col, xshort ? x.value.size : y.value.size);
                return res != 0 ? res < 0 : xshort;
              });
  }
};

// Key is a C string, sort alphabetically case insensitive.
struct ustr_type {
  const char* col;
  size_t      size;
};
template<>
struct header_traits<ustr_type> {
  typedef header_type<ustr_type> type;

  static type name(const char* start, const char* end, const char* col) {
    ustr_type s = { col, strcspn(std::min(col, end), " \n\t") };
    return type(s, start, end);
  }

  static void sort(std::vector<type>& headers) {
    std::sort(headers.begin(), headers.end(),
              [](const type& x, const type& y) -> bool {
                const bool xshort = x.value.size < y.value.size;
                const int res = strncasecmp(x.value.col, y.value.col, xshort ? x.value.size : y.value.size);
                return res != 0 ? res < 0 : xshort;
              });
  }
};

// No key. Sort randomly
struct random_type { };
template<>
struct header_traits<random_type> {
  typedef header_type<random_type> type;

  static type name(const char* start, const char* end, const char* col) {
    return type(random_type(), start, end);
  }

  static void sort(std::vector<type>& headers) {
    std::random_shuffle(headers.begin(), headers.end());
  }
};

// Return a pointer to the nb-th space separated token in str. str is
// not modified. Returns NULL if less than nb columns.
const char* find_token(uint32_t nb, const char* str) {
  static const char* space = " \t\n";
  for(uint32_t i = 1; i < nb && *str && *str != '\n'; ++i) {
    str += strcspn(str, space);
    str += strspn(str, space);
  }
  return (*str && *str != '\n') ? str : nullptr;
}

template<typename T>
void parse_headers(const char* const start, const char* const end, std::vector<header_type<T>>& headers) {
  const char* current = start;
  // Ignore up to first header
  if(*current != '>') {
    for(current = strchr(current, '>'); current; current = strchr(current + 1, '>'))
      if(*(current - 1) == '\n') break;
  }

  // At the beginning of the loop, current points to the start of the
  // header '>'
  for(const char* next = current + 1; current; current = next) {
    for(next = strchr(next + 1, '>'); next; next = strchr(next + 1, '>'))
      if(*(next - 1) == '\n') break;
    const char* key = find_token(args.key_arg, current + 1);
    const size_t eol = strcspn(key, "\n");
    headers.push_back(header_traits<T>::name(current, next ? next : end - 1, key + std::min(eol, (size_t)args.character_arg - 1)));
  }
}

template <typename T>
static int sort_mmap(const sort_cmdline& args) {
  //  static_assert(std::is_nothrow_move_constructible<close_fd>::value, "Close_fd nothrow move");
  std::vector<close_fd>       fds;
  std::vector<header_type<T>> headers;

  for(auto path : args.file_arg) {
    fds.push_back(open_mmap(path));
    auto& current = fds.back();
    parse_headers<T>(current.ptr, current.ptr + current.size, headers);
  }

  header_traits<T>::sort(headers);

  for(const auto& header : headers)
    std::cout.write(header.start, header.end - header.start);

  return EXIT_SUCCESS;
}

int sort_main(int argc, char *argv[]) {
  args.parse(argc, argv);
  if(args.file_arg.empty()) {
    args.file_arg.push_back("/dev/stdin");
    if(isatty(0))
      std::cerr << "Warning: reading from terminal" << std::endl;
  }

  if(args.numeric_sort_flag)
    return sort_mmap<int64_t>(args);
  else if(args.random_sort_flag)
    return sort_mmap<random_type>(args);
  else if(args.ignore_case_flag)
    return sort_mmap<ustr_type>(args);
  else
    return sort_mmap<str_type>(args);
}
