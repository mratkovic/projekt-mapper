/*
 * lcskpp.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 */

#include <algorithm>

#include <metrics_algorithm/lcskpp.h>
#include "util/util_structures.h"
#include "util/fenwick.h"

using namespace std;

uint32_t LCSkpp::calcLCSkpp(uint32_t k,
                            vector<pair<uint32_t, uint32_t> > &result,
                            vector<pair<uint32_t, uint32_t> > &elements) {

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
  Fenwick<pair<uint32_t, uint32_t> > maxColDp(n);

  vector<uint32_t> dp(elements.size());
  vector<int> recon(elements.size());
  vector<int> continues(elements.size(), -1);

  // find pairs continuing each other
  if (k > 1) {
    vector<pair<uint32_t, uint32_t> >::iterator it;
    vector<pair<uint32_t, uint32_t> >::iterator prevIt;

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

  reconstructLCSkpp(elements, k, recon, bestIndex, lcskppLen, result);
  return lcskppLen;
}

uint32_t LCSkpp::calcLCSkppSlow(
    uint32_t k, vector<pair<uint32_t, uint32_t>> &reconstruction,
    vector<pair<uint32_t, uint32_t>> &matches) {

  if (matches.empty()) {
    reconstruction.clear();
    return 0;
  }

  int n = matches.size();
  vector<int> dp(n);
  vector<int> recon(n);

  int bestIndex = 0;
  int lcskLen = 0;

  for (int i = 0; i < n; ++i) {
    dp[i] = k;
    recon[i] = -1;

    int rowEnd = matches[i].first + k - 1;
    int colEnd = matches[i].second + k - 1;
    int primDiagI = rowEnd - colEnd;
    int secDiagI = (rowEnd + colEnd) / 2;

    for (int j = i - 1; j >= 0; --j) {
      if (matches[j].first + k <= matches[i].first
          && matches[j].second + k <= matches[i].second) {
        // 1) Uzimam cijeli match interval i nastavljam neki
        // match koji je ranije vec zavrsio.
        if (dp[j] + (int) k > dp[i]) {
          dp[i] = dp[j] + k;
          recon[i] = j;
        }
      } else {
        // 2) Nastavak po istoj dijagonali.
        int rowEndJ = matches[j].first + k - 1;
        int colEndJ = matches[j].second + k - 1;
        int primDiagJ = rowEndJ - colEndJ;
        int secDiagJ = (rowEndJ + colEndJ) / 2;

        int extend = secDiagI - secDiagJ;
        if (primDiagI == primDiagJ && secDiagI > secDiagJ && extend < (int) k) {
          if (dp[j] + extend > dp[i]) {
            dp[i] = dp[j] + extend;
            recon[i] = j;
          }
        }

      }
    }

    if (dp[i] > lcskLen) {
      bestIndex = i;
      lcskLen = dp[i];
    }
  }

  reconstructLCSkpp(matches, k, recon, bestIndex, lcskLen, reconstruction);
  return lcskLen;
}

void LCSkpp::reconstructLCSkpp(
    vector<pair<uint32_t, uint32_t>> &elements, uint32_t k,
    vector<int> &prevIndex, int lastIndex, int lcskLen,
    vector<pair<uint32_t, uint32_t>> &reconstruction) {

  reconstruction.clear();
  reconstruction.reserve(lcskLen);

  int index = lastIndex;
  while (index != -1) {
    int refEndIndex = elements[index].first + k - 1;
    int readEndIndex = elements[index].second + k - 1;

    int prev = prevIndex[index];
    uint32_t howManyElements;

    bool takeWhole = prev == -1;
    if (prev != -1 && elements[prev].first + k <= elements[index].first
        && elements[prev].second + k <= elements[index].second) {
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

