/*
 * suffix_array.cpp
 *
 *  Created on: Dec 27, 2014
 *      Author: marko
 */

#include "suffix_array.h"

SuffixArray::SuffixArray(const char* text, size_t textLen)
    : text_(text),
      textLen_(textLen) {

  suffix_array_.clear();
  suffix_array_.resize(textLen);
  divsufsort((const sauchar_t *) text, &suffix_array_[0], textLen);
}
SuffixArray::SuffixArray(FILE* in, const char* text, size_t textLen) {
  uint32_t len;
  suffix_array_.clear();
  fread(&len, sizeof(len), 1, in);
  printf("%d -- %d\n", len, textLen);
  assert(len == textLen);

  text_ = text;
  textLen_ = textLen;
  suffix_array_.resize(textLen_);
  fread(&suffix_array_[0], sizeof(std::vector<saidx_t>::value_type),
        suffix_array_.size(), in);
}
SuffixArray::~SuffixArray() {
  text_ = NULL;
}

const int *SuffixArray::search(const char* pattern, int length,
                               int* numOfSolution) {
  int firstIndex;
  *numOfSolution = sa_search((const sauchar_t *) text_, textLen_,
                             (const sauchar_t *) pattern, length,
                             &suffix_array_[0], textLen_, &firstIndex);

  if (*numOfSolution == -1) {
    return numOfSolution;
  } else {
    return &suffix_array_[firstIndex];
  }
}

//const int* SuffixArray::iterativeSearch(const char* pattern, int length,
//                                        int startLen, int* numOfSolutions,
//                                        const int solUpperLimit,
//                                        const int solLowerLimit,
//                                        int* finalLen) {
//  int firstIndex = -1;
//  int prevFirstIndex;
//
//  int prevNumOfSolutions;
//  *numOfSolutions = -1;
//
//  int currentLen = startLen - 1;
//
//  // TODO: +-+-
//  do {
//    currentLen++;
//    prevNumOfSolutions = *numOfSolutions;
//    prevFirstIndex = firstIndex;
//    *numOfSolutions = sa_search((const sauchar_t *) text_, textLen_,
//                                (const sauchar_t *) pattern, currentLen,
//                                &suffix_array_[0], textLen_, &firstIndex);
//
//  } while (*numOfSolutions > solUpperLimit && currentLen <= length);
//
//  *finalLen = currentLen;
//  if (*numOfSolutions < solLowerLimit && prevNumOfSolutions != -1) {
//    *numOfSolutions = prevNumOfSolutions;
//    firstIndex = prevFirstIndex;
//    *finalLen = currentLen - 1;
//  }
//
//  if (*numOfSolutions == -1) {
//    return numOfSolutions;
//  } else {
//    return &suffix_array_[firstIndex];
//  }
//
//}

const int* SuffixArray::iterativeSearch(const char* pattern, int length,
                                        int startLen, int* numOfSolutions,
                                        const int solUpperLimit,
                                        const int solLowerLimit,
                                        int* finalLen) {
  int firstIndex = -1;
  int prevFirstIndex;

  int prevNumOfSolutions;
  *numOfSolutions = -1;

  int currentLen = startLen - 1;

  // TODO: +-+-  exp pa lin
  do {
    currentLen++;
    prevNumOfSolutions = *numOfSolutions;
    prevFirstIndex = firstIndex;
    *numOfSolutions = sa_search((const sauchar_t *) text_, textLen_,
                                (const sauchar_t *) pattern, currentLen,
                                &suffix_array_[0], textLen_, &firstIndex);

  } while (*numOfSolutions > solUpperLimit && currentLen <= length);

  *finalLen = currentLen;
  if (*numOfSolutions < solLowerLimit && prevNumOfSolutions != -1) {
    *numOfSolutions = prevNumOfSolutions;
    firstIndex = prevFirstIndex;
    *finalLen = currentLen - 1;
  }

  if (*numOfSolutions == -1) {
    return numOfSolutions;
  } else {
    return &suffix_array_[firstIndex];
  }

}

void SuffixArray::saveSuffixArray(FILE* out) {
  fwrite(&textLen_, sizeof(textLen_), 1, out);
  fwrite(&suffix_array_[0], sizeof(std::vector<saidx_t>::value_type),
         suffix_array_.size(), out);
}

uint32_t SuffixArray::size() {
  return textLen_;
}
const char* SuffixArray::text() {
  return text_;
}
