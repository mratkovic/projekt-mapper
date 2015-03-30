/*
 * lcsk.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 */

#include "lcsk.h"
#include <algorithm>
#include <cstdio>
#include <map>

#include "../util/fenwick.h"

using namespace std;

struct Event {
  uint32_t first;
  uint32_t second;
  bool isStart;  // start or end;
  uint32_t index;
  Event(uint32_t first, uint32_t second, bool isStart, uint32_t index)
      : first(first),
        second(second),
        isStart(isStart),
        index(index) {
  }

};

bool operator<(const Event &a, const Event &b) {
  if (a.first != b.first)
    return a.first < b.first;

  if (a.second != b.second)
    return a.second < b.second;

  return a.isStart < b.isStart;
}

struct Event3 {
  uint32_t first;
  uint32_t second;
  uint32_t k;
  bool isStart;  // start or end;
  uint32_t index;
  Event3(uint32_t first, uint32_t second, uint32_t k, bool isStart,
         uint32_t index)
      : first(first),
        second(second),
        k(k),
        isStart(isStart),
        index(index) {
  }

};

bool operator<(const Event3 &a, const Event3 &b) {
  if (a.first != b.first)
    return a.first < b.first;

  if (a.second != b.second)
    return a.second < b.second;

  return a.isStart < b.isStart;
}
bool operator==(const Event3 &a, const Event3 &b) {

  return a.first == b.first && a.second == b.second && (a.isStart == b.isStart);
}

void reconstructLCSk(vector<pair<uint32_t, uint32_t> >& elements, uint32_t k,
                     vector<int>& prevIndex, int lastIndex, int lcskLen,
                     vector<pair<uint32_t, uint32_t> >& reconstruction) {

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

void reconstructLCSkpp(vector<pair<uint32_t, uint32_t> >& elements, uint32_t k,
                       vector<int>& prevIndex, int lastIndex, int lcskLen,
                       vector<pair<uint32_t, uint32_t> >& reconstruction) {

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

uint32_t LCSk::calcLCSk(uint32_t k, vector<pair<uint32_t, uint32_t> >& result,
                        vector<pair<uint32_t, uint32_t> >& elements) {

  if (elements.empty()) {
    result.clear();
    return 0;
  }

  vector<Event> events;
  events.reserve(2 * elements.size());
  uint32_t n = 0;

  for (uint32_t i = 0; i < elements.size(); ++i) {
    pair<uint32_t, uint32_t> element = elements[i];

    events.push_back(Event(element.first, element.second, true, i));
    events.push_back(Event(element.first + k, element.second + k, false, i));

    n = max(n, element.second + k);
  }
  sort(events.begin(), events.end());

  // Indexed by column, first:dp value, second:index in elements
  Fenwick<pair<uint32_t, uint32_t> > maxColDp(n);

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

uint32_t LCSk::calcLCSkpp(uint32_t k, vector<pair<uint32_t, uint32_t> >& result,
                          vector<pair<uint32_t, uint32_t> >& elements) {

  if (elements.empty()) {
    result.clear();
    return 0;
  }

  vector<Event> events;
  events.reserve(2 * elements.size());
  uint32_t n = 0;

  for (uint32_t i = 0; i < elements.size(); ++i) {
    pair<uint32_t, uint32_t> element = elements[i];

    events.push_back(Event(element.first, element.second, true, i));
    events.push_back(Event(element.first + k, element.second + k, false, i));

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

  for (vector<Event>::iterator event = events.begin(); event != events.end();
      ++event) {
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

uint32_t reconstructLCSpp(vector<Triplet>& elements, vector<int>& prevIndex,
                          int lastIndex, int lcskLen,
                          vector<pair<uint32_t, uint32_t> > &reconstruction) {

  int score = 0;

  reconstruction.clear();
  reconstruction.reserve(lcskLen);

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
//
//      int dx = elements[index].first - elements[prev].first - prevK - 1;
//      int dy = elements[index].second - elements[prev].second - prevK - 1;
//
//      score += N * min(dx, dy);
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

    //score += howManyElements * E;
    for (uint32_t j = 0; j < howManyElements; ++j) {
      reconstruction.push_back(make_pair(refEndIndex - j, readEndIndex - j));
    }
    index = prevIndex[index];
  }

  reverse(reconstruction.begin(), reconstruction.end());

  return std::max<int>(0, score);
}

uint32_t LCSk::estimateBeginingPosFromLCSk(
    vector<pair<uint32_t, uint32_t> >& reconstruction) {
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

uint32_t LCSk::calcLCSkppSlow(
    uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >& reconstruction,
    std::vector<std::pair<uint32_t, uint32_t> >& matches) {

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
        if (dp[j] + k > dp[i]) {
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
        if (primDiagI == primDiagJ && secDiagI > secDiagJ && extend < k) {
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

uint32_t LCSk::getScoreFromRecon(
    uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >& reconstruction) {
  int score = 0;

  if (reconstruction.size() == 0) {
    return 0;
  }
  for (int i = 0; i < reconstruction.size() - 1; ++i) {
    int dx = reconstruction[i + 1].first - reconstruction[i].first - 1;
    int dy = reconstruction[i + 1].second - reconstruction[i].second - 1;

    score += N * (dx + dy) + E;

  }
  score += E;

  uint32_t sc = std::max<int>(0, score);
  return sc;
}

uint32_t LCSk::calcLCSkppPenalty(uint32_t k,
                                 vector<pair<uint32_t, uint32_t> >& result,
                                 vector<pair<uint32_t, uint32_t> >& elements) {

  if (elements.empty()) {
    result.clear();
    return 0;
  }

  vector<Event> events;
  events.reserve(2 * elements.size());
  uint32_t n = 0;

  for (uint32_t i = 0; i < elements.size(); ++i) {
    pair<uint32_t, uint32_t> element = elements[i];

    events.push_back(Event(element.first, element.second, true, i));
    events.push_back(Event(element.first + k, element.second + k, false, i));

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

  for (vector<Event>::iterator event = events.begin(); event != events.end();
      ++event) {
    int index = event->index;

    if (event->isStart) {
      pair<int, int> max = maxColDp.getMax(event->second);
      dp[index] = k * E;
      recon[index] = -1;

      if (max.first > 0) {
        dp[index] = max.first + N * event->first + N * event->second + E * k;
        recon[index] = max.second;
      }

    } else {
//      if (continues[index] != -1) {
//        if (dp[continues[index]] + 1 > dp[index]) {
//          dp[index] = dp[continues[index]] + E;
//          recon[index] = continues[index];
//        }
//      }
      maxColDp.updateMax(
          event->second,
          make_pair(dp[index] - N * (event->first + event->second + 2 * k),
                    index));

      if (dp[index] > lcskppLen) {
        lcskppLen = dp[index];
        bestIndex = index;
      }
    }
  }

  reconstructLCSkpp(elements, k, recon, bestIndex, lcskppLen, result);
  return lcskppLen;
}

uint32_t LCSk::calcLCSpp(std::vector<std::pair<uint32_t, uint32_t> >& result,
                         std::vector<Triplet>& elements) {

  if (elements.empty()) {
    result.clear();
    return 0;
  }

  vector<Event3> events;
  events.reserve(2 * elements.size());
  uint32_t n = 0;

  for (uint32_t i = 0; i < elements.size(); ++i) {
    Triplet element = elements[i];

    uint32_t k = element.third;
    events.push_back(Event3(element.first, element.second, k, true, i));
    events.push_back(
        Event3(element.first + k, element.second + k, k, false, i));

    n = max(n, element.second + k);
  }
  sort(events.begin(), events.end());

  // Indexed by column, first:dp value, second:index in elements
  Fenwick<pair<uint32_t, uint32_t> > maxColDp(n);

  vector<uint32_t> dp(elements.size());
  vector<int> recon(elements.size());
  vector<int> continues(elements.size(), -1);

  // find pairs continuing each other

  vector<Triplet>::iterator it;
  vector<Triplet>::iterator prevIt;

  for (it = elements.begin(); it != elements.end(); ++it) {
    Triplet prevElement = Triplet(it->first - 1, it->second - 1, 0);
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
  return lcskppLen;
}

