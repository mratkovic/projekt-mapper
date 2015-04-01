/*
 * util_structures.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_METRICS_ALGORITHM_UTIL_STRUCTURES_H_
#define SRC_METRICS_ALGORITHM_UTIL_STRUCTURES_H_

struct event_t {
  uint32_t first;
  uint32_t second;
  bool isStart;
  uint32_t index;

  event_t(uint32_t first, uint32_t second, bool isStart, uint32_t index)
      : first(first),
        second(second),
        isStart(isStart),
        index(index) {
  }

  bool operator<(const event_t &b) const {
    if (first != b.first)
      return first < b.first;

    if (second != b.second)
      return second < b.second;

    return isStart < b.isStart;
  }

};

struct eventK_t {
  uint32_t first;
  uint32_t second;
  bool isStart;
  uint32_t index;

  uint32_t k;

  eventK_t(uint32_t first, uint32_t second, uint32_t k, bool isStart,
           uint32_t index)
      : first(first),
        second(second),
        k(k),
        isStart(isStart),
        index(index) {
  }

};

template<typename T>
struct triplet_t {
  T first;
  T second;
  T third;

  triplet_t(T first, T second, T third)
      : first(first),
        second(second),
        third(third) {
  }

};

inline bool operator<(const triplet_t<uint32_t> &a,
                      const triplet_t<uint32_t> &b) {
  if (a.first != b.first) {
    return a.first < b.first;
  }

  return a.second < b.second;

}

#endif /* SRC_METRICS_ALGORITHM_UTIL_STRUCTURES_H_ */
