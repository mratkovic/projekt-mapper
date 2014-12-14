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
#include "Mapping.h"

#include "BioUtil.h"

#define KMER_K 20
#define LIS true
#define WINDOW_SIZE 2

int estimateBeginingPosFromLIS(std::vector<std::pair<int, int> >& positions, std::vector<int>& lis) {
	std::map<int, int> possiblePosCntr;

	for (size_t i = 0; i < lis.size(); ++i) {
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
		int kmerStart, bool complement) {
	char array[KMER_K];
	for (int i = 0; i < KMER_K; ++i) {
		array[i] = read->getData(kmerStart + i, complement);
	}
	sa->findStartingPositions(array, KMER_K, kmerStart, positions);
}
void mapReadToSuffixArray(bioutil::Read* read, SuffixArray* sa) {
	std::vector<std::pair<int, int> > pos;
	mapReadToSuffixArray(read, sa, false);
	// komplement???
	mapReadToSuffixArray(read, sa, true);
}

void mapReadToSuffixArray(bioutil::Read* read, SuffixArray* sa, bool complement) {
	std::vector<std::pair<int, int> > pos;
	for (int i = 0; i < read->getDataLen() - KMER_K; ++i) {
		getKhmerPositions(read, sa, pos, i, complement);
	}
	std::sort(pos.begin(), pos.end());
//	printf("%d\n", pos.size());
//	for (int i = 0; i < pos.size(); ++i) {
//		printf(" %d     %d ---- %d\n", i, pos[i].first, pos[i].second);
//	}

	if (pos.size() == 0) {
		return;
	}
	// sweep
	int startIndex = 0, endIndex = 0, lastPosition = -1000000;
	while (true) {
		for (; startIndex < pos.size() && pos[startIndex].first < lastPosition + read->getDataLen(); ++startIndex) {
			// pomakni pocetak prozora za duljinu reada od proslog pocetka
			;
		}

		if (startIndex < pos.size()) {
			lastPosition = pos[startIndex].first;
		} else {
			break;
		}
		int windowEnd = pos[startIndex].first + WINDOW_SIZE * read->getDataLen();
		for (; endIndex < pos.size() && pos[endIndex].first < windowEnd; ++endIndex) {
			// prozor dux max WINDOW SIZE * duzina reada
			;
		}
		assert(startIndex <= endIndex);

		if (LIS) {
			runLIS(startIndex, endIndex, pos, read, complement);
		} else {
			runKLCS(startIndex, endIndex, pos, read, complement);
		}
	}
}

void runLIS(int startIndex, int endIndex, std::vector<std::pair<int, int> > &pos, bioutil::Read* read,
		bool complemented) {
	std::vector<std::pair<int, int> > lisData;
	for (int i = startIndex; i < endIndex; ++i) {
		lisData.push_back(pos[i]);
	}
	std::vector<int> lisResult;
	calcLIS(&lisResult, lisData);
	int beginPos = estimateBeginingPosFromLIS(pos, lisResult);
	bioutil::Mapping mapping(lisResult.size(), beginPos, beginPos + read->getDataLen(), complemented);
	read->addMapping(mapping);
}
void runKLCS(int startIndex, int endIndex, std::vector<std::pair<int, int> > &pos, bioutil::Read* read,
		bool complemented) {
	// zuzicev algoritam iz pdf-a
	// spora implementacija
	std::vector<std::pair<int, int> > matchPairs(pos.begin() + startIndex, pos.begin() + endIndex);

	if (matchPairs.empty()) {
		return;
	}
	int startInSequence = matchPairs[0].first;
	int maxSeq = 0, maxRead = 0;

	// translatirat na 0
	for (size_t i = 0; i < matchPairs.size(); ++i) {
		matchPairs[i].first -= startInSequence;
		maxSeq = std::max(maxSeq, matchPairs[i].first);
		maxRead = std::max(maxRead, matchPairs[i].second);
	}

	int klcsResult = 0;
	int n = matchPairs.size();
	std::vector<int> dp(n);
	std::vector<int> reconstruct(n);
	int bestFind = 0;

	for (int i = 0; i < n; ++i) {
		dp[i] = KMER_K;
		reconstruct[i] = -1;

		int rowEnd = matchPairs[i].first + KMER_K - 1;
		int colEnd = matchPairs[i].second + KMER_K - 1;

		int primDiagI = rowEnd - colEnd;
		int secDiagI = (rowEnd + colEnd) / 2;

		for (int j = i - 1; j >= 0; --j) {
			if (matchPairs[j].first + KMER_K <= matchPairs[i].first
					&& matchPairs[j].second + KMER_K <= matchPairs[i].second) {

				// jedan chunk od k i drugi chunk od k
				// nema spajanja
				if (dp[j] + KMER_K > dp[i]) {
					dp[i] = dp[j] + KMER_K;
					reconstruct[i] = j;
				}
			} else {
				int rowEndJ = matchPairs[j].first + KMER_K - 1;
				int colEndJ = matchPairs[j].second + KMER_K - 1;

				int primDiagJ = rowEndJ - colEndJ;
				int secDiagJ = (rowEndJ + colEndJ) / 2;

				if (primDiagI == primDiagJ && secDiagI > secDiagJ && secDiagI - secDiagJ < KMER_K) {
					int addition = secDiagI - secDiagJ;
					if (dp[j] + addition > dp[i]) {
						dp[i] = dp[j] + addition;
						reconstruct[i] = j;
					}
				}
			}
		}

		if (dp[i] > klcsResult) {
			bestFind = i;
			klcsResult = dp[i];
		}
	}

	std::vector<std::pair<int, int> > reconstruction;
	reconstruction.reserve(klcsResult);
}

int calcEditDistanceNaive(char* pattern1, char* pattern2, int len1, int len2) {
	std::vector<int> vi(len2 + 1, 0);
	std::vector<std::vector<int> > dp;
	dp.insert(dp.begin(), len1 + 1, vi);

	for (int i = 0; i <= len1; ++i) {
		dp[i][0] = i;
	}

	for (int j = 0; j <= len2; ++j) {
		dp[0][j] = j;
	}

	for (int i = 1; i <= len1; ++i) {
		for (int j = 1; j <= len2; ++j) {
			int match = dp[i - 1][j - 1] + (pattern1[i - 1] == pattern2[j - 1] ? 0 : 1);
			int insert = dp[i][j - 1] + 1;
			int del = dp[i - 1][j] + 1;
			dp[i][j] = std::min(match, insert);
			dp[i][j] = std::min(dp[i][j], del);
		}
	}
	return dp[len1][len2];
}

int calcEditDistanceNaiveMem(char* pattern1, char* pattern2, int len1, int len2) {
	// len1, pattern1 kraci
	if (len1 > len2) {
		int pom;
		char* pomPtr;

		pom = len1;
		len1 = len2;
		len2 = pom;

		pomPtr = pattern1;
		pattern1 = pattern2;
		pattern2 = pomPtr;
	}

	int pom = 0;
	int* dp = new int[len1 + 1];
	for (int i = 0; i <= len1; ++i) {
		dp[i] = i;
	}

	for (int i = 1; i <= len1; ++i) {
		for (int j = 1; j <= len2; ++j) {
			if (j == 1) {
				pom = dp[j - 1];
				dp[j - 1] = i;
			}
			int match = pom + (pattern1[i - 1] == pattern2[j - 1] ? 0 : 1);
			int insert = dp[j - 1] + 1;
			int del = dp[j] + 1;

			pom = dp[j];

			dp[j] = std::min(match, insert);
			dp[j] = std::min(dp[j], del);
		}
	}
	pom = dp[len2];
	delete[] dp;
	return pom;
}
