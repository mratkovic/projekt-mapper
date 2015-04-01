/*
 * lcskpp.h
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 */
#ifndef SRC_METRICS_ALGORITHM_LCSKPP_H_
#define SRC_METRICS_ALGORITHM_LCSKPP_H_

#include <vector>
#include <stdint.h>

class LCSkpp {
 public:

  /**

   * @param k parameter k of LCSk++
   * @result  vector for keeping reconstruction of LCSk++. It contains pairs of integers that
   *          correspond with positions in strings.
   * @elements vector filled with pairs (i, j) that represent starting positions of common substrings
   *            of length k. i is starting position in first string, j starting position in second string
   * @return  value of LCSk++ between two strings
   */
  static uint32_t calcLCSkpp(
      uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >& result,
      std::vector<std::pair<uint32_t, uint32_t> >& elements);

  static uint32_t estimateBeginingPosFromLCSkpp(
      std::vector<std::pair<uint32_t, uint32_t> >& reconstruction);

  static uint32_t calcLCSkppSlow(
      uint32_t k, std::vector<std::pair<uint32_t, uint32_t>> &reconstruction,
      std::vector<std::pair<uint32_t, uint32_t>> &matches);

  static void reconstructLCSkpp(
      std::vector<std::pair<uint32_t, uint32_t>> &elements, uint32_t k,
      std::vector<int> &prevIndex, int lastIndex, int lcskLen,
      std::vector<std::pair<uint32_t, uint32_t>> &reconstruction);
};

#endif /* SRC_METRICS_ALGORITHM_LCSKPP_H_ */
