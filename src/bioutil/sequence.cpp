/*
 * sequence.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include <cassert>
#include <vector>
#include "sequence.h"
#include "../external/kseq.h"
#include "../util/utility_functions.h"
#include "bioutil.h"

namespace bioutil {

KSEQ_INIT(int, read)

Sequence::Sequence() {
  //data_ = info_ = comment_ = NULL;
  data_ = NULL;
  dataLen_ = 0;
  basesAsInt = false;
  numOfSequences_ = 0;
}
Sequence::~Sequence() {
  clear();
}

void Sequence::clear() {
  if (data_) {
    delete[] data_;
    data_ = 0;
  }

  for (auto it = info_.begin(); it != info_.end(); ++it) {
    delete *it;
  }
  info_.clear();

  dataLen_ = 0;
}

const char* Sequence::data() {
  return data_;
}
const char* Sequence::info(uint32_t index) {
  return info_[index];
}

uint32_t Sequence::dataLen() {
  return dataLen_;
}

uint32_t Sequence::numOfSequences() {
  return numOfSequences_;
}

void Sequence::readSequencesFromFASTA(FILE* fastaIn) {
  clear();

  fpos_t begin;
  fgetpos(fastaIn, &begin);
  kseq_t *seq = kseq_init(fileno(fastaIn));

  dataLen_ = 0;
  numOfSequences_ = 0;

  while (kseq_read(seq) >= 0) {
    uint32_t size = seq->name.l + 1;
    char* info = new char[size];
    memcpy(info, seq->name.s, seq->name.l);
    info[seq->name.l] = 0;
    info_.push_back(info);

    uint32_t dataLen = seq->seq.l;
    char* data = new char[dataLen];
    memcpy(data, seq->seq.s, dataLen);

    ++numOfSequences_;

    dataLen_ += dataLen;
    seqEndIndex_.push_back(dataLen_);

  }

  assert(seqEndIndex_.size() == info_.size());
  assert(seqEndIndex_.size() == numOfSequences_);
  data_ = new char[dataLen_];

  kseq_destroy(seq);

  fsetpos(fastaIn, &begin);
  seq = kseq_init(fileno(fastaIn));

  uint32_t currentSize = 0;
  while (kseq_read(seq) >= 0) {
    uint32_t dataLen = seq->seq.l;
    memcpy(data_ + currentSize, seq->seq.s, dataLen);

    currentSize += dataLen;
  }

}

void Sequence::readSingleSequenceFromFASTA(FILE* fastaIn) {
  clear();
  kseq_t *seq = kseq_init(fileno(fastaIn));
  assert(kseq_read(seq) >= 0);

  numOfSequences_ = 1;
  uint32_t size = seq->name.l + 1;
  info_.push_back(new char[size]);
  memcpy(info_[0], seq->name.s, seq->name.l);
  info_[0][seq->name.l] = 0;

  dataLen_ = seq->seq.l;
  data_ = new char[dataLen_];
  memcpy(data_, seq->seq.s, dataLen_);

  seqEndIndex_.push_back(dataLen_);

  kseq_destroy(seq);
}
void Sequence::allBasesToSmallInt() {
  if (basesAsInt) {
    return;
  }

  for (uint32_t i = 0; i < dataLen_; ++i) {
    data_[i] = baseToInt(data_[i]);
  }
  basesAsInt = true;
}
void Sequence::allBasesToLetters() {
  if (!basesAsInt) {
    return;
  }

  for (uint32_t i = 0; i < dataLen_; ++i) {
    data_[i] = intToBase(data_[i]);
  }
  basesAsInt = false;
}

uint32_t Sequence::sequenceIndex(uint32_t positionGlobal) {
  // binary search

  uint32_t lo, hi, mid;
  lo = 0;
  hi = numOfSequences_ - 1;

  while (lo < hi) {
    mid = lo + (hi - lo) / 2;
    uint32_t lastPos = seqEndIndex_[mid];

    if (positionGlobal < lastPos) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }

  return lo;
}

uint32_t Sequence::seqLen(uint32_t index) {
  if (index == 0) {
    return seqEndIndex_[0];
  } else {
    return seqEndIndex_[index] - seqEndIndex_[index - 1];
  }
}

uint32_t Sequence::positionInSeq(uint32_t positionGlobal) {
  uint32_t seqIndex = sequenceIndex(positionGlobal);

  if (seqIndex == 0) {
    return positionGlobal;
  } else {
    return positionGlobal - seqEndIndex_[seqIndex - 1];
  }
  return 0;
}

}

