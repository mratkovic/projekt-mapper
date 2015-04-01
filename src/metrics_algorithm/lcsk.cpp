/*
 * lcsk.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 */

#include <algorithm>
#include <cstdio>
#include <map>

#include "lcsk.h"
#include "util/fenwick.h"
#include "util/util_structures.h"

using namespace std;

uint32_t LCSk::calcLCSk(uint32_t k, vector<pair<uint32_t, uint32_t>>& result,
                        vector<pair<uint32_t, uint32_t>>& elements) {

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

  uint32_t lcskLen = 0;
  uint32_t bestIndex = 0;

  for (auto event = events.begin(); event != events.end(); ++event) {
    int index = event->index;

    if (event->isStart) {
      pair<int, int> max = maxColDp.getMax(event->second);
      dp[index] = k;
      recon[index] = -1;

      if (max.first > 0) {
        dp[index] = max.first + k;
        recon[index] = max.second;
      }

    } else {

      maxColDp.updateMax(event->second, make_pair(dp[index], index));

      if (dp[index] > lcskLen) {
        lcskLen = dp[index];
        bestIndex = index;
      }
    }
  }

  reconstructLCSk(elements, k, recon, bestIndex, lcskLen, result);
  return lcskLen;
}

uint32_t LCSk::estimateBeginingPosFromLCSk(
    vector<pair<uint32_t, uint32_t>>& reconstruction) {
  map<uint32_t, uint32_t> possiblePosCntr;

  for (auto match : reconstruction) {
    uint32_t realPos = match.first;
    uint32_t kmerStartingPos = match.second;
    ++possiblePosCntr[max<uint32_t>(0, realPos - kmerStartingPos)];
  }

  uint32_t maxVotes = 0;
  uint32_t begin = reconstruction[0].first - reconstruction[0].second;
  map<int, int>::iterator it;

  for (auto pair : possiblePosCntr) {
    if (pair.second > maxVotes) {
      maxVotes = pair.second;
      begin = pair.first;
    }
  }
  return begin;
}

void LCSk::reconstructLCSk(vector<pair<uint32_t, uint32_t>> &elements,
                           uint32_t k, vector<int>& prevIndex, int lastIndex,
                           int lcskLen,
                           vector<pair<uint32_t, uint32_t>> &reconstruction) {

  reconstruction.clear();
  reconstruction.reserve(lcskLen);

  int index = lastIndex;
  while (index != -1) {
    int refEndIndex = elements[index].first + k - 1;
    int readEndIndex = elements[index].second + k - 1;

    for (uint32_t j = 0; j < k; ++j) {
      reconstruction.push_back(make_pair(refEndIndex - j, readEndIndex - j));
    }

    index = prevIndex[index];
  }

  reverse(reconstruction.begin(), reconstruction.end());
}


