/*
 * suffix_array.h
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

 public:
  SuffixArray(const char* text, size_t textLen);
  SuffixArray(FILE *in, const char* text, size_t textLen);
  ~SuffixArray();

  /**
   * @param pattern to search for
   * @param length length of pattern
   * @param numOfSolutions number of found patterns in the text
   * @returns the pointer to the index of the first occurrence
   */
  const int* search(const char* pattern, int length, int* numOfSolutions);
  void saveSuffixArray(FILE *out);

  uint32_t size();
  const char* text();

 private:
  const char* text_;
  uint32_t textLen_;
  std::vector<saidx_t> suffix_array_;
};

#endif

