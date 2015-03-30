/*
 * triplet.h
 *
 *  Created on: Mar 30, 2015
 *      Author: marko
 */

#ifndef SRC_UTIL_TRIPLET_H_
#define SRC_UTIL_TRIPLET_H_

struct Triplet {
  uint32_t first;
  uint32_t second;
  uint32_t third;

  Triplet(uint32_t first, uint32_t second, uint32_t third)
      : first(first),
        second(second),
        third(third) {
  }

};

inline bool operator<(const Triplet &a, const Triplet &b) {
  if (a.first != b.first)
    return a.first < b.first;

  return a.second < b.second;

}

#endif /* SRC_UTIL_TRIPLET_H_ */
