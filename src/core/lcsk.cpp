/*
 * lcsk.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 */

#include "lcsk.h"
#include <algorithm>
#include <cstdio>

#include "../util/fenwick.h"

struct Event {
	uint32_t first;
	uint32_t second;
	bool isStart; // start or end;
	uint32_t index;
	Event(uint32_t first, uint32_t second, bool isStart, uint32_t index) :
			first(first), second(second), isStart(isStart), index(index) {
	}

};

bool operator<(const Event &a, const Event &b) {
	if (a.first != b.first)
		return a.first < b.first;

	if (a.second != b.second)
		return a.second < b.second;

	return a.isStart < b.isStart;
}

void reconstructLCSk(std::vector<std::pair<uint32_t, uint32_t> >& elements, uint32_t k, std::vector<int>& prevIndex,
		int lastIndex, int lcskLen, std::vector<std::pair<uint32_t, uint32_t> >* reconstruction) {

	reconstruction->clear();
	reconstruction->reserve(lcskLen);

	int index = lastIndex;
	while (index != -1) {
		int refEndIndex = elements[index].first + k - 1;
		int readEndIndex = elements[index].second + k - 1;

		for (uint32_t j = 0; j < k; ++j) {
			reconstruction->push_back(std::make_pair(refEndIndex - j, readEndIndex - j));
		}

		index = prevIndex[index];
	}

	std::reverse(reconstruction->begin(), reconstruction->end());
}

uint32_t LCSk::calcLCSk(uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >* result,
		std::vector<std::pair<uint32_t, uint32_t> >& elements) {

	if (elements.empty()) {
		result->clear();
		return 0;
	}

	std::vector<Event> events;
	events.reserve(2 * elements.size());
	uint32_t n = 0;

	for (uint32_t i = 0; i < elements.size(); ++i) {
		std::pair<uint32_t, uint32_t> element = elements[i];

		events.push_back(Event(element.first, element.second, true, i));
		events.push_back(Event(element.first + k, element.second + k, false, i));

		n = std::max(n, element.second + k);
	}
	sort(events.begin(), events.end());

	// Indexed by column, first:dp value, second:index in elements
	Fenwick<std::pair<uint32_t, uint32_t> > maxColDp(n);

	std::vector<uint32_t> dp(elements.size());
	std::vector<int> recon(elements.size());

	uint32_t lcskLen = 0;
	uint32_t bestIndex = 0;

	for (auto event = events.begin(); event != events.end(); ++event) {
		int index = event->index;

		if (event->isStart) {
			std::pair<int, int> max = maxColDp.getMax(event->second);
			dp[index] = k;
			recon[index] = -1;

			if (max.first > 0) {
				dp[index] = max.first + k;
				recon[index] = max.second;
			}

		} else {

			maxColDp.updateMax(event->second, std::make_pair(dp[index], index));

			if (dp[index] > lcskLen) {
				lcskLen = dp[index];
				bestIndex = index;
			}
		}
	}

	reconstructLCSk(elements, k, recon, bestIndex, lcskLen, result);

	return lcskLen;
}

uint32_t LCSk::calcLCSkpp(uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >* result,
		std::vector<std::pair<uint32_t, uint32_t> >& elements) {

	if (elements.empty()) {
		result->clear();
		return 0;
	}

	std::vector<Event> events;
	events.reserve(2 * elements.size());
	uint32_t n = 0;

	for (uint32_t i = 0; i < elements.size(); ++i) {
		std::pair<uint32_t, uint32_t> element = elements[i];

		events.push_back(Event(element.first, element.second, true, i));
		events.push_back(Event(element.first + k, element.second + k, false, i));

		n = std::max(n, element.second + k);
	}
	sort(events.begin(), events.end());

	// Indexed by column, first:dp value, second:index in elements
	Fenwick<std::pair<uint32_t, uint32_t> > maxColDp(n);

	std::vector<uint32_t> dp(elements.size());
	std::vector<int> recon(elements.size());
	std::vector<int> continues(elements.size(), -1);

	// find pairs continuing each other
	if (k > 1) {
		std::vector<std::pair<uint32_t, uint32_t> >::iterator it;
		std::vector<std::pair<uint32_t, uint32_t> >::iterator prevIt;

		for (it = elements.begin(); it != elements.end(); ++it) {
			std::pair<uint32_t, uint32_t> prevElement = std::make_pair(it->first - 1, it->second - 1);
			prevIt = std::lower_bound(elements.begin(), elements.end(), prevElement);
			if (*prevIt == prevElement) {
				continues[it - elements.begin()] = prevIt - elements.begin();
			}
		}
	}

	uint32_t lcskppLen = 0;
	uint32_t bestIndex = 0;

	for (std::vector<Event>::iterator event = events.begin(); event != events.end(); ++event) {
		int index = event->index;

		if (event->isStart) {
			std::pair<int, int> max = maxColDp.getMax(event->second);
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
			maxColDp.updateMax(event->second, std::make_pair(dp[index], index));

			if (dp[index] > lcskppLen) {
				lcskppLen = dp[index];
				bestIndex = index;
			}
		}
	}

	reconstructLCSk(elements, k, recon, bestIndex, lcskppLen, result);

	return lcskppLen;
}

