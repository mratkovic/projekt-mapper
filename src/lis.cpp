/*
 * lis.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 */

#include "lis.h"
#include "fenwick.h"

#include <map>
#include <algorithm>

namespace bioutil {

int Lis::estimateBeginingPosFromLIS(std::vector<std::pair<uint32_t, uint32_t> >& positions, std::vector<int>& lis) {
	std::map<int, int> possiblePosCntr;

	for (uint32_t i = 0; i < lis.size(); ++i) {
		int relativePos = lis[i];
		int realPos = positions[relativePos].first;
		int kmerStartingPos = positions[relativePos].second;

		++possiblePosCntr[std::max(0, realPos - kmerStartingPos)];
	}

	int maxVotes = -1;
	int begin = -1;
	std::map<int, int>::iterator it;
	for (it = possiblePosCntr.begin(); it != possiblePosCntr.end(); ++it) {
		if (it->second > maxVotes) {
			maxVotes = it->second;
			begin = it->first;
		}
	}
	return begin;
}

void Lis::calcLIS(std::vector<int>* result, const std::vector<std::pair<uint32_t, uint32_t> >& elements) {
	uint32_t n = elements.size();
	if (n == 0) {
		return;
	}
	int* dp = new int[n];
	int* dpPath = new int[n];

	for (uint32_t i = 0; i < n; ++i) {
		dp[i] = 0;
		dpPath[i] = -1;
	}

	uint32_t maxValue = elements[0].second;
	for (uint32_t i = 0; i < elements.size(); ++i) {
		maxValue = std::max(maxValue, elements[i].second);
	}

	Fenwick<std::pair<int, int> > fnwk(maxValue);

	for (uint32_t i = 0; i < n; ++i) {
		std::pair<int, int> max = fnwk.getMax(elements[i].second - 1);
		if (max.first) {
			dp[i] = max.first + 1;
			dpPath[i] = max.second;
		} else {
			dp[i] = 1;
		}
		fnwk.updateMax(elements[i].second, std::make_pair(dp[i], i));
	}

	int max = -1;
	int pos = -1;
	for (uint32_t i = 0; i < n; ++i) {
		if (max < dp[i]) {
			pos = i;
			max = dp[i];
		}
	}
	reconstructLIS(result, pos, dpPath, n);

	delete[] dp;
	delete[] dpPath;
}
void Lis::reconstructLIS(std::vector<int>* result, int lastPos, int* reconstructionTable, int n) {
	result->clear();
	result->reserve(n);

	while (lastPos != -1) {
		result->push_back(lastPos);
		lastPos = reconstructionTable[lastPos];
	}

	reverse(result->begin(), result->end());
}

} /* namespace bioutil */
