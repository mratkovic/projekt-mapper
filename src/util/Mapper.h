/*
 * Mapper.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#ifndef SRC_UTIL_MAPPER_H_
#define SRC_UTIL_MAPPER_H_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <dirent.h>
#include <unistd.h>
#include <utility>
#include <vector>

#include "SuffixArray.h"
#include "Read.h"
#include "Fenwick.h"

void calcLIS(std::vector<int>* result, const std::vector<std::pair<int, int> >& elements);
void reconstructLIS(std::vector<int>* result, int lastPos, int* dpPath, int n);
int estimateBeginingPosFromLIS(std::vector<std::pair<int, int> >& positions, std::vector<int>& lis);

void mapReadToSuffixArray(bioutil::Read* read, SuffixArray* sa);
void getKhmerPositions(bioutil::Read* read, SuffixArray* sa, std::vector<std::pair<int, int> > &positions,
		int kmerStart);
void mapReadToSuffixArray(bioutil::Read* read, SuffixArray* sa, bool complement);
void getPositions(bioutil::Read* read, SuffixArray* sa, bool complement);

void runLIS(int startIndex, int endIndex, std::vector<std::pair<int, int> > &pos, bioutil::Read* read, bool complement);
void runKLCS(int startIndex, int endIndex, std::vector<std::pair<int, int> > &pos, bioutil::Read* read,
		bool complement);

int calcEditDistanceNaive(char* pattern1, char* pattern2, int len1, int len2);
int calcEditDistanceNaiveMem(char* pattern1, char* pattern2, int len1, int len2);
#endif
