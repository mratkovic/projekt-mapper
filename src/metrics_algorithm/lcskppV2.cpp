/*
 * lcskppV2.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#include <algorithm>

#include "lcskppV2.h"
#include "util/fenwick.h"
#include <cassert>

using namespace std;

bool operator<(const eventK_t &a, const eventK_t &b) {
  if (a.first != b.first)
    return a.first < b.first;

  if (a.second != b.second)
    return a.second < b.second;

  if (a.index != b.index)
    return a.index < b.index;

  if (a.isStart != b.isStart)
    return a.isStart < b.isStart;

  // LOLOLOL
  if (a.k != b.k)
    return a.k < b.k;

}

bool operator==(const eventK_t &a, const eventK_t &b) {

  return a.first == b.first && a.second == b.second && (a.isStart == b.isStart);
}

uint32_t LCSkppV2::calcLCSkpp(
    std::vector<std::pair<uint32_t, uint32_t>> &result,
    std::vector<triplet_t<uint32_t>> &elements) {

  if (elements.empty()) {
    result.clear();
    return 0;
  }

  vector<eventK_t> events;
  events.reserve(2 * elements.size());
  uint32_t n = 0;

  for (uint32_t i = 0; i < elements.size(); ++i) {
    triplet_t<uint32_t> element = elements[i];

    uint32_t k = element.third;
    events.push_back(eventK_t(element.first, element.second, k, true, i));
    events.push_back(
        eventK_t(element.first + k, element.second + k, k, false, i));

    n = max(n, element.second + k);
  }
  sort(events.begin(), events.end());

  // Indexed by column, first:dp value, second:index in elements
  Fenwick<pair<uint32_t, uint32_t>> maxColDp(n);

  vector<uint32_t> dp(elements.size());
  vector<int> recon(elements.size());
  vector<int> continues(elements.size(), -1);

  // find pairs continuing each other
  vector<triplet_t<uint32_t>>::iterator it;
  vector<triplet_t<uint32_t>>::iterator prevIt;

  for (it = elements.begin(); it != elements.end(); ++it) {
    triplet_t<uint32_t> prevElement = triplet_t<uint32_t>(it->first - 1,
                                                          it->second - 1, 0);
    prevIt = lower_bound(elements.begin(), elements.end(), prevElement);
    if (prevIt->first == prevElement.first
        && prevIt->second == prevElement.second) {
      continues[it - elements.begin()] = prevIt - elements.begin();
    }
  }

  uint32_t lcskppLen = 0;
  uint32_t bestIndex = 0;

  for (auto event = events.begin(); event != events.end(); ++event) {
    int index = event->index;

    if (event->isStart) {
      pair<int, int> max = maxColDp.getMax(event->second);
      dp[index] = event->k;
      recon[index] = -1;

      if (max.first > 0) {
        dp[index] = max.first + event->k;
        recon[index] = max.second;
      }

    } else {

      if (continues[index] != -1) {
        if (dp[continues[index]] + 1 > dp[index]) {
          dp[index] = dp[continues[index]] + 1;
          recon[index] = continues[index];
        }
      }
      maxColDp.updateMax(event->second, make_pair(dp[index], index));

      if (dp[index] > lcskppLen) {
        lcskppLen = dp[index];
        bestIndex = index;
      }
    }
  }

  reconstructLCSpp(elements, recon, bestIndex, lcskppLen, result);
  result.erase(unique(result.begin(), result.end()), result.end());
  //assert(result.size() == lcskppLen);
  return result.size();
}

void LCSkppV2::reconstructLCSpp(
    vector<triplet_t<uint32_t>> &elements, vector<int> &prevIndex,
    int lastIndex, int lcskLen,
    vector<pair<uint32_t, uint32_t>> &reconstruction) {

  reconstruction.clear();
  //reconstruction.reserve(lcskLen);

  int index = lastIndex;
  while (index != -1) {

    int k = elements[index].third;
    int refEndIndex = elements[index].first + k - 1;
    int readEndIndex = elements[index].second + k - 1;

    int prev = prevIndex[index];
    int prevK = elements[prev].third;

    uint32_t howManyElements;

    bool takeWhole = prev == -1;
    if (prev != -1 && elements[prev].first + prevK <= elements[index].first
        && elements[prev].second + prevK <= elements[index].second) {
      takeWhole = true;
    }

    if (takeWhole) {
      howManyElements = k;
    } else {
      int curr_secondary_diag = (elements[index].first + elements[index].second)
          / 2;
      int prev_secondary_diag = (elements[prev].first + elements[prev].second)
          / 2;
      howManyElements = curr_secondary_diag - prev_secondary_diag;
    }

    for (uint32_t j = 0; j < howManyElements; ++j) {
      reconstruction.push_back(make_pair(refEndIndex - j, readEndIndex - j));
    }
    index = prevIndex[index];
  }

  reverse(reconstruction.begin(), reconstruction.end());
}

