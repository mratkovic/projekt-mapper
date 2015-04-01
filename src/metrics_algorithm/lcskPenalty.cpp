/*
 * lcskPenalty.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#include <metrics_algorithm/lcskPenalty.h>
#include "util/fenwick.h"
#include "lcskpp.h"
#include "util/util_structures.h"

#include <algorithm>

using namespace std;

uint32_t LCSkPenalty::getScoreFromRecon(
    uint32_t k, vector<pair<uint32_t, uint32_t>> &reconstruction,
    int match_score, int manhattan_penalty) {

  int score = 0;
  if (reconstruction.size() == 0) {
    return 0;
  }
  for (uint32_t i = 0; i < reconstruction.size() - 1; ++i) {
    int dx = reconstruction[i + 1].first - reconstruction[i].first - 1;
    int dy = reconstruction[i + 1].second - reconstruction[i].second - 1;

    score += manhattan_penalty * (dx + dy) + match_score;

  }
  score += match_score;

  return max<int>(0, score);
}

uint32_t LCSkPenalty::calcLCSkppPenalty(
    uint32_t k, vector<pair<uint32_t, uint32_t>> &result,
    vector<pair<uint32_t, uint32_t>> &elements, int match_score,
    int manhattan_penalty) {

  if (elements.empty()) {
    result.clear();
    return 0;
  }

  vector<event_t> events;
  events.reserve(2 * elements.size());
  uint32_t n = 0;

  for (uint32_t i = 0; i < elements.size(); ++i) {
    pair<uint32_t, uint32_t> element = elements[i];

    events.push_back(event_t(element.first, element.second, true, i));
    events.push_back(event_t(element.first + k, element.second + k, false, i));

    n = max(n, element.second + k);
  }
  sort(events.begin(), events.end());

// Indexed by column, first:dp value, second:index in elements
  Fenwick<pair<uint32_t, uint32_t>> maxColDp(n);

  vector<uint32_t> dp(elements.size());
  vector<int> recon(elements.size());
  vector<int> continues(elements.size(), -1);

// find pairs continuing each other
  if (k > 1) {
    vector<pair<uint32_t, uint32_t>>::iterator it;
    vector<pair<uint32_t, uint32_t>>::iterator prevIt;

    for (it = elements.begin(); it != elements.end(); ++it) {
      pair<uint32_t, uint32_t> prevElement = make_pair(it->first - 1,
                                                       it->second - 1);
      prevIt = lower_bound(elements.begin(), elements.end(), prevElement);
      if (*prevIt == prevElement) {
        continues[it - elements.begin()] = prevIt - elements.begin();
      }
    }
  }

  uint32_t lcskppLen = 0;
  uint32_t bestIndex = 0;

  for (vector<event_t>::iterator event = events.begin(); event != events.end();
      ++event) {
    int index = event->index;

    if (event->isStart) {
      pair<int, int> max = maxColDp.getMax(event->second);
      dp[index] = k * match_score;
      recon[index] = -1;

      if (max.first > 0) {
        dp[index] = max.first + manhattan_penalty * event->first
            + manhattan_penalty * event->second + match_score * k;
        recon[index] = max.second;
      }

    } else {
//      if (continues[index] != -1) {
//        if (dp[continues[index]] + 1 > dp[index]) {
//          dp[index] = dp[continues[index]] + E;
//          recon[index] = continues[index];
//        }
//      }

      int value = dp[index]
          - manhattan_penalty * (event->first + event->second + 2 * k);

      maxColDp.updateMax(event->second, make_pair(value, index));

      if (dp[index] > lcskppLen) {
        lcskppLen = dp[index];
        bestIndex = index;
      }
    }
  }

  LCSkpp::reconstructLCSkpp(elements, k, recon, bestIndex, lcskppLen, result);
  return lcskppLen;
}

