#include "util.h"

#include <chrono>

using namespace std;

namespace util {

// stack of times of start_timer calls
thread_local static vector<chrono::time_point<std::chrono::steady_clock>> start_times;

void start_timer()
{
    start_times.push_back(chrono::steady_clock::now());
}

double stop_timer()
{
    auto stop_time = chrono::steady_clock::now();
    double diff_nano = chrono::duration_cast<chrono::nanoseconds>(stop_time - start_times.back()).count();
    start_times.pop_back();
    return diff_nano / 1.e9;
}

Summary Summary::operator*(double x) const
{
    return { min * x, max * x, avg * x  };
}

}

namespace std {

static ListFormat current_list_format = ListFormat::plain;

void set_list_format(ListFormat format)
{
    current_list_format = format;
}

ListFormat get_list_format()
{
    return current_list_format;
}

std::ostream& operator<<(std::ostream& os, util::Summary s)
{
    return os << s.min << " - " << s.max << " (avg " << s.avg << ")";
}

}
