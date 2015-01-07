/*
 * Sequence.h
 *
 *  Created on: Dec 27, 2014
 *      Author: marko
 */

#ifndef DTRA_SUFFIX_ARRAY_DATABASE
#define DTRA_SUFFIX_ARRAY_DATABASE

#include <divsufsort.h>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <vector>
#include <stdint.h>

class SuffixArray {

private:
	const char* text_;
	uint32_t textLen_;
	std::vector<saidx_t> suffix_array_;

public:
	SuffixArray(const char* text, size_t textLen);
	SuffixArray(FILE *in, const char* text, size_t textLen);
	~SuffixArray();

	/**
	 * @param pattern to search for
	 * @param patternLen
	 * @param numSolutions number of found patterns in the text
	 * @returns the pointer to the index of the first occurrence; the second occurrence is on *(pointer+1) and so on
	 */
	const int* search(const char* pattern, int patternLen, int* numOfSolutions);
	void saveSuffixArray(FILE *out);

	uint32_t size();
	const char* text();
};

#endif

