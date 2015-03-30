/*
 * lcsk.h
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 *  Static class that defines method for calculation of LCSk between two strings.
 *  LCSk is number of common substrings of length k between two given strings.
 */

#ifndef SRC_LCSK_H_
#define SRC_LCSK_H_

#include <vector>
#include <stdint.h>
#include "../util/triplet.h"

#define E 10
#define N -1

class LCSk {

 public:
  /**
   * Method that is used for calculation of LCSk between two given strings with already
   * filled vector of positions of every common substring of length k in those two strings.
   *
   * @param k parameter k of LCSk
   * @result vector for keeping reconstruction of LCSk. It contains pairs of integers that
   * correspond with positions in strings.
   * @elements vector filled with pairs (i, j) that represent starting positions of common substrings
   * of length k. i is starting position in first string, j starting position in second string
   * @return value of LCSk between two strings
   */
  static uint32_t calcLCSk(
      uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >& result,
      std::vector<std::pair<uint32_t, uint32_t> >& elements);

  /**
   * Method that is used for calculation of LCSk++ between two given strings with already
   * filled vector of positions of every common substring of length k in those two strings.
   *
   * @param k parameter k of LCSk++
   * @result vector for keeping reconstruction of LCSk++. It contains pairs of integers that
   * correspond with positions in strings.
   * @elements vector filled with pairs (i, j) that represent starting positions of common substrings
   * of length k. i is starting position in first string, j starting position in second string
   * @return value of LCSk++ between two strings
   */
  static uint32_t calcLCSkpp(
      uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >& result,
      std::vector<std::pair<uint32_t, uint32_t> >& elements);

  static uint32_t calcLCSpp(std::vector<std::pair<uint32_t, uint32_t> >& result,
                            std::vector<Triplet>& elements);

  // nepomaze za velik postotak greske
  static uint32_t estimateBeginingPosFromLCSk(
      std::vector<std::pair<uint32_t, uint32_t> >& reconstruction);

  static uint32_t calcLCSkppSlow(
      uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >& result,
      std::vector<std::pair<uint32_t, uint32_t> >& elements);

  static uint32_t calcLCSkppPenalty(
      uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >& result,
      std::vector<std::pair<uint32_t, uint32_t> >& elements);

  static uint32_t getScoreFromRecon(
      uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >& reconstruction);
};

#endif /* SRC_LCSK_H_ */
