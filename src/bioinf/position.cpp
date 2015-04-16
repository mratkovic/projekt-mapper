/*
 * mapping.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#include <bioinf/position.h>
#include <cstring>

Position::Position() {
  score_ = secondary_score_ = start_ = end_ = 0;
  complement_ = false;
  cigar_ = 0;
}

Position::Position(uint32_t score, uint32_t secondaryScore, uint32_t start,
                   uint32_t end, bool complemented, const char* cigarStr,
                   uint32_t cigarLen)
    : score_(score),
      secondary_score_(secondaryScore),
      start_(start),
      end_(end),
      complement_(complemented),
      cigar_(0) {

      cigar(cigarStr, cigarLen);
}

bool Position::isComplement() {
  return complement_;
}
void Position::setComplement(bool complemented) {
  complement_ = complemented;
}

uint32_t Position::score() const {
  return score_;
}
uint32_t Position::secondaryScore() const {
  return secondary_score_;
}

uint32_t Position::end() {
  return end_;
}

uint32_t Position::start() const {
  return start_;
}

void Position::end(uint32_t end) {
  end_ = end;
}
void Position::start(uint32_t start) {
  start_ = start;
}

void Position::score(uint32_t score) {
  score_ = score;
}

const char* Position::cigar() {
  return cigar_;
}
void Position::cigar(const char* cigar, uint32_t len) {
  if (cigar == NULL) {
    return;
  }

  cigar_ = new char[len + 1];
  memcpy(cigar_, cigar, len);
  cigar_[len] = 0;
}

Position::~Position() {
  if (cigar_) {
    delete[] cigar_;
    cigar_ = 0;
  }
}
