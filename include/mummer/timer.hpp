#ifndef __MUMMER_TIMER_H__
#define __MUMMER_TIMER_H__


#ifdef MEASURE_TIME
#include <chrono>
#include <atomic>
#include <ostream>
#include <iomanip>


namespace timer {

typedef std::chrono::steady_clock::time_point time_point;
typedef std::chrono::microseconds             duration;
inline time_point time() { return std::chrono::steady_clock::now(); }
inline duration since(const time_point& t) {
  return std::chrono::duration_cast<duration>(time() - t);
}

struct add_time_to {
  duration&  m_timer;
  time_point m_start;

  add_time_to(duration& t) : m_timer(t), m_start(time()) { }
  ~add_time_to() { m_timer += since(m_start); }
};

struct global {
  std::atomic<uint64_t> m_time;
  global() : m_time(0) { }
  global& operator+=(const duration& rhs) {
    m_time += rhs.count();
    return *this;
  }
  operator duration() const { return duration(m_time.load()); }
};

inline std::ostream& operator<<(std::ostream& os, duration x) {
  auto h = std::chrono::duration_cast<std::chrono::hours>(x);
  x -= h;
  auto m = std::chrono::duration_cast<std::chrono::minutes>(x);
  x -= m;
  auto s = std::chrono::duration_cast<std::chrono::seconds>(x);
  x -= s;
  auto i = std::chrono::duration_cast<std::chrono::milliseconds>(x);
  return os << h.count() << ':'
            << std::setfill('0') << std::setw(2) << m.count() << ':'
            << std::setfill('0') << std::setw(2) << s.count() << '.'
            << std::setfill('0') << std::setw(3) << i.count();
}

struct scope {
  time_point    m_start;
  std::ostream& m_os;
  const char*   m_msg;
  scope(std::ostream& os = std::clog) : m_start(time()), m_os(os), m_msg(nullptr) { }
  scope(const char* msg, std::ostream& os = std::clog) : m_start(time()), m_os(os), m_msg(msg) { }
  ~scope() {
    if(m_msg)
      m_os << m_msg << ' ';
    m_os << since(m_start) << std::endl;
  }
};
} // namespace timer
#define TIME_SCOPE(msg) timer::scope scope_timer ## __COUNTER__(msg);
#define TIME_FUNCTION timer::scope function_timer ## __COUNTER__(__PRETTY_FUNCTION__);
#else
#include <ostream>
namespace timer {
struct timer { };
struct add_time_to {
  add_time_to(const timer& t) { }
};
struct global {
  template<typename T>
  global& operator+=(const T& rhs) { return *this; }
  std::string to_hms() const { return std::string(); }
};
inline std::ostream& operator<<(std::ostream& os, const global& rhs) {
  return os;
}

struct scope {
  scope(std::ostream& os = std::clog) { }
  scope(const char* msg, std::ostream& os = std::clog) { }
};
} // namespace timer
#define TIME_SCOPE(msg)
#define TIME_FUNCTION

#endif

#endif /* __MUMMER_TIMER_H__ */
