/*
 * lcskPenalty.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_METRICS_ALGORITHM_LCSKPENALTY_H_
#define SRC_METRICS_ALGORITHM_LCSKPENALTY_H_

#include <vector>
#include <cstdint>
#include <utility>

class LCSkPenalty {
 public:
  static uint32_t getScoreFromRecon(
      uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >& reconstruction,
      int match_score, int manhattan_penalty);

  uint32_t calcLCSkppPenalty(
      uint32_t k, std::vector<std::pair<uint32_t, uint32_t>> &result,
      std::vector<std::pair<uint32_t, uint32_t>> &elements, int match_score,
      int manhattan_penalty);

};

#endif /* SRC_METRICS_ALGORITHM_LCSKPENALTY_H_ */
