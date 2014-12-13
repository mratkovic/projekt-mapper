/*
 * Mapper.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#include <cassert>
#include <cstdio>
#include <algorithm>
#include <utility>
#include <vector>
#include <map>
#include "Mapper.h"

#include "BioUtil.h"

#define KMER_K 20
#define LIS true
#define WINDOW_SIZE 2

int estimateBeginingPosFromLIS(std::vector<std::pair<int, int> >& positions, std::vector<int>& lis) {
	std::map<int, int> possiblePos;

	for (size_t i = 0; i < lis.size(); ++i) {
		int relativePos = lis[i];
		int realPos = positions[relativePos].first;
		int kmerStartingPos = positions[relativePos].second;

		++possiblePos[std::max(0, realPos - kmerStartingPos)];
	}

	int maxVotes = -1;
	int begin = -1;
	std::map<int, int>::iterator it;
	for (it = possiblePos.begin(); it != possiblePos.end(); ++it) {
		if (it->second > maxVotes) {
			maxVotes = it->second;
			begin = it->first;
		}
	}
	printf("Procijena pozicije - %d\n", begin);
	return begin;
}

void calcLIS(std::vector<int>* result, const std::vector<std::pair<int, int> >& elements) {
	int n = elements.size();
	if (n == 0) {
		return;
	}
	int* dp = new int[n];
	int* dpPath = new int[n];
	for (int i = 0; i < n; ++i) {
		dp[i] = 0;
		dpPath[i] = -1;
	}

	int maxValue = -1;
	for (int i = 0; i < elements.size(); ++i) {
		maxValue = std::max(maxValue, elements[i].second);
	}

	Fenwick<std::pair<int, int> > fnwk(maxValue);

	for (size_t i = 0; i < n; ++i) {
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
	for (int i = 0; i < n; ++i) {
		if (max < dp[i]) {
			pos = i;
			max = dp[i];
		}
	}
	reconstructLIS(result, pos, dpPath, n);
}
void reconstructLIS(std::vector<int>* result, int lastPos, int* reconstructionTable, int n) {
	result->clear();
	result->reserve(n);
	while (lastPos != -1) {
		result->push_back(lastPos);
		lastPos = reconstructionTable[lastPos];
	}
	reverse(result->begin(), result->end());
}

void getKhmerPositions(bioutil::Read* read, SuffixArray* sa, std::vector<std::pair<int, int> > &positions,
		int kmerStart) {
	char intArray[KMER_K];
	for(int i = 0; i < KMER_K; ++i) {
		intArray[i] = baseToInt(read->getData(kmerStart+i));
	}
	sa->findStartingPositions(intArray, KMER_K, kmerStart, positions);
}
int getPositionInGeneFromSuffixArray(bioutil::Read* read, SuffixArray* sa) {
	std::vector<std::pair<int, int> > pos;
	printf("iden pozicije trazit\n");
	for (int i = 0; i < read->getDataLen() - KMER_K; ++i) {
		getKhmerPositions(read, sa, pos, i);
	}
	std::sort(pos.begin(), pos.end());
	// sweep
	int startIndex = 0, endIndex = 0, lastPosition = -1;
	while (true) {
		while (startIndex < pos.size() && lastPosition + read->getDataLen() > pos[startIndex].first) {
			++startIndex;
		}

		if (startIndex < pos.size()) {
			lastPosition = pos[startIndex].first;
		} else {
			break;
		}
		int windowEnd = pos[startIndex].first + WINDOW_SIZE * read->getDataLen();
		for (; endIndex < pos.size() && pos[endIndex].first < windowEnd; ++endIndex)
			;
		assert(startIndex <= endIndex);
		if (LIS) {
			printf("LISSS\n");
			std::vector<int> lisResult;
			calcLIS(&lisResult, pos);
			int beginPos = estimateBeginingPosFromLIS(pos, lisResult);

			printf("BEGINNNNN %d\n", beginPos);
			char info[250];
			sprintf(info, "%ul, %d, %ul", lisResult.size(), beginPos, beginPos + read->getDataLen());
			read->setMappingInfo(info);
			return beginPos;
		} else {
			std::vector<int> result;
			calcLIS(&result, pos);
		}
	}

	return -1;
}

