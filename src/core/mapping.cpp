/*
 * mapping.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#include <cstring>
#include "mapping.h"

Mapping::Mapping() {
  score_ = start_ = end_ = 0;
  complement_ = false;
  cigar_ = 0;
}
Mapping::Mapping(uint32_t score, uint32_t start, uint32_t end, bool complement,
                 const char* cigarStr, uint32_t cigarLen)
    : score_(score),
      start_(start),
      end_(end),
      complement_(complement),
      cigar_(0) {
  cigar(cigarStr, cigarLen);
}

bool Mapping::isComplement() {
  return complement_;
}
void Mapping::setComplement(bool complemented) {
  complement_ = complemented;
}

uint32_t Mapping::score() const {
  return score_;
}

uint32_t Mapping::end() {
  return end_;
}

uint32_t Mapping::start() {
  return start_;
}

void Mapping::end(uint32_t end) {
  end_ = end;
}
void Mapping::start(uint32_t start) {
  start_ = start;
}

const char* Mapping::cigar() {
  return cigar_;
}
void Mapping::cigar(const char* cigar, uint32_t len) {
  if (cigar == NULL) {
    return;
  }

  cigar_ = new char[len + 1];
  memcpy(cigar_, cigar, len);
  cigar_[len] = 0;
}

Mapping::~Mapping() {
  if (cigar_) {
    delete[] cigar_;
    cigar_ = 0;
  }
}
