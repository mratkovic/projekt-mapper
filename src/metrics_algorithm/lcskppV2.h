/*
 * lcskppV2.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_METRICS_ALGORITHM_LCSKPPV2_H_
#define SRC_METRICS_ALGORITHM_LCSKPPV2_H_

#include <vector>
#include <utility>
#include <cstdint>
#include "util/util_structures.h"

class LCSkppV2 {
 public:

  static uint32_t calcLCSkpp(std::vector<std::pair<uint32_t, uint32_t>> &result,
                             std::vector<triplet_t<uint32_t>> &elements);

  static uint32_t estimateBeginingPosFromLCSkpp(
      std::vector<std::pair<uint32_t, uint32_t> >& reconstruction);

 private:
  static void reconstructLCSpp(
      std::vector<triplet_t<uint32_t>> &elements, std::vector<int> &prevIndex,
      int lastIndex, int lcskLen,
      std::vector<std::pair<uint32_t, uint32_t>> &reconstruction);

};
#endif /* SRC_METRICS_ALGORITHM_LCSKPPV2_H_ */
