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
