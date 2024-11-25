#pragma once

#include <vector>
#include <algorithm>
#include <ostream>

namespace util {

// start new time measurement
void start_timer();
// returns time in seconds since last unconsumed start_timer call and consumes it
double stop_timer();

// sort vector and remove duplicate elements
template<typename T>
void make_set(std::vector<T> &v)
{
    size_t v_size = v.size();
    if (v_size == 0)
        return;
    std::sort(v.begin(), v.end());
    size_t last_distinct = 0;
    for (size_t next = 1; next < v_size; next++)
        if (v[next] != v[last_distinct])
        {
            last_distinct++;
            std::swap(v[next], v[last_distinct]);
        }
    v.resize(last_distinct + 1);
}

// remove elements in set from v
// set must be sorted
template<typename T>
void remove_set(std::vector<T> &v, const std::vector<T> set)
{
    assert(is_sorted(set.cbegin(), set.cend()));
    if (v.empty() || set.empty())
        return;
    std::erase_if(v, [&set](T value) { return std::binary_search(set.cbegin(), set.cend(), value); });
}

struct Summary
{
    double min;
    double max;
    double avg;
    Summary operator*(double x) const;
};

template<typename T, class Map>
Summary summarize(const std::vector<T> &v, Map f)
{
    Summary summary = {};
    if (!v.empty())
        summary.min = f(v[0]);
    for (const T& e : v)
    {
        double x = f(e);
        summary.avg += x;
        if (x < summary.min)
            summary.min = x;
        if (x > summary.max)
            summary.max = x;
    }
    if (!v.empty())
        summary.avg /= v.size();
    return summary;
}

// compute total number of elements in vector of collections
template<typename T>
size_t size_sum(const std::vector<T> &v)
{
    size_t sum = 0;
    for (const T &x : v)
        sum += x.size();
    return sum;
}

// extract size values of vector of collection
template<typename T>
std::vector<size_t> sizes(const std::vector<T> &v)
{
    std::vector<size_t> s;
    for (const T &x : v)
        s.push_back(x.size());
    return s;
}

template<typename T>
T random(const std::vector<T> &v)
{
    assert(v.size() > 0);
    return v[rand() % v.size()];
}

template<typename T>
class min_bucket_queue
{
    std::vector<std::vector<T>> buckets;
    size_t min_bucket; // minimum non-empty bucket
public:
    min_bucket_queue() : min_bucket(0) {}
    void push(T value, size_t bucket)
    {
        if (empty() || min_bucket > bucket)
            min_bucket = bucket;
        if (buckets.size() <= bucket)
            buckets.resize(bucket + 1);
        buckets[bucket].push_back(value);
    }
    bool empty() const
    {
        return min_bucket >= buckets.size();
    }
    T pop() {
        assert(min_bucket < buckets.size() && !buckets[min_bucket].empty());
        T top = buckets[min_bucket].back();
        buckets[min_bucket].pop_back();
        // skip empty buckets
        while (min_bucket < buckets.size() && buckets[min_bucket].empty())
            min_bucket++;
        return top;
    }
};

} // util

namespace std {

enum class ListFormat { plain, indexed };

void set_list_format(ListFormat format);
ListFormat get_list_format();

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &v)
{
    if (v.empty())
        return os << "[]";
    if (get_list_format() == ListFormat::indexed)
    {
        os << "[0:" << v[0];
        for (size_t i = 1; i < v.size(); i++)
            os << ',' << i << ":" << v[i];
    }
    else
    {
        os << "[" << v[0];
        for (size_t i = 1; i < v.size(); i++)
            os << ',' << v[i];
    }
    return os << ']';
}

template <typename A, typename B>
std::ostream& operator<<(std::ostream& os, const std::pair<A,B> &p)
{
    return os << "(" << p.first << "," << p.second << ")";
}

std::ostream& operator<<(std::ostream& os, util::Summary s);

} // std
