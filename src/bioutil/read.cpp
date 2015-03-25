/*
 * read.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include "read.h"
#include "../util/utility_functions.h"

#define LINE_SIZE 30000
#define KEEP_FACTOR 1l

using namespace bioutil;
Read::Read() {
  id_ = data_ = optional_identifier_ = quality_ = 0;
  dataLen_ = 0;
  basesAsInt_ = false;
  //mapping_ = new Mapping;
}
Read::~Read() {
  clear();
}

void Read::addMapping(double score, uint32_t start, uint32_t end,
                      bool isComplement, const char* cigar, uint32_t cigarLen) {
  Mapping* m = new Mapping(score, start, end, isComplement, cigar, cigarLen);
  mappings_.insert(m);

  while (mappings_.size() > 1
      && ((*mappings_.rbegin())->score() / (*mappings_.begin())->score()
          > KEEP_FACTOR)) {
    delete *mappings_.begin();
    mappings_.erase(mappings_.begin());
  }
}

Mapping* Read::bestMapping(uint32_t index) {
  std::multiset<Mapping*, ptr_compare<Mapping> >::reverse_iterator it =
      mappings_.rbegin();
  uint32_t cntr = 0;

  for (; it != mappings_.rend(); ++it, ++cntr) {
    if (cntr == index) {
      return *it;
    }
  }
  return NULL;
}
std::multiset<Mapping*, ptr_compare<Mapping> >& Read::mappings() {
  return mappings_;
}
uint32_t Read::mappingsSize() {
  return mappings_.size();
}

const char* Read::data() {
  return data_;
}
const char* Read::id() {
  return id_;
}

uint32_t Read::dataLen() {
  return dataLen_;
}

void Read::clear() {
  if (data_) {
    delete[] data_;
    data_ = 0;
  }
  if (id_) {
    delete[] id_;
    id_ = 0;
  }
  if (optional_identifier_) {
    delete[] optional_identifier_;
    optional_identifier_ = 0;
  }
  if (quality_) {
    delete[] quality_;
    quality_ = 0;
  }

  std::multiset<Mapping*, ptr_compare<Mapping> >::iterator it =
      mappings_.begin();
  uint32_t cntr = 0;

  for (; it != mappings_.end(); ++it, ++cntr) {
    delete *it;
  }
  dataLen_ = 0;

}

bool Read::readNextFromFASTQ(kseq_t* seq) {
  if (kseq_read(seq) < 0) {
    return false;
  }

  uint32_t size = seq->name.l + 1;
  id_ = new char[size];
  memcpy(id_, seq->name.s, seq->name.l);
  id_[size - 1] = 0;

  dataLen_ = seq->seq.l;
  data_ = new char[dataLen_ + 1];
  memcpy(data_, seq->seq.s, dataLen_);
  data_[dataLen_] = 0;

  quality_ = new char[dataLen_ + 1];
  memcpy(quality_, seq->qual.s, dataLen_);
  quality_[dataLen_] = 0;

  return true;
}

void Read::allBasesToSmallInt() {
  if (basesAsInt_) {
    return;
  }

  for (uint32_t i = 0; i < dataLen_; ++i) {
    data_[i] = baseToInt(data_[i]);
  }
  basesAsInt_ = true;
}
void Read::allBasesToLetters() {
  if (!basesAsInt_) {
    return;
  }
  for (uint32_t i = 0; i < dataLen_; ++i) {
    data_[i] = intToBase(data_[i]);
  }
  basesAsInt_ = false;
}

Read* Read::getReverseComplement() {
  Read *rev = new Read();
  rev->dataLen_ = dataLen_;
  rev->data_ = new char[dataLen_ + 1];

  for (uint32_t i = 0; i < dataLen_; ++i) {
    if (basesAsInt_) {
      rev->data_[i] = getACGTComplementAsSmallInt(data_[dataLen_ - 1 - i]);
    } else {
      rev->data_[i] = getACGTComplement(data_[dataLen_ - 1 - i]);
    }
  }
  rev->data_[dataLen_] = 0;
  return rev;

}
void Read::printReadSAM(FILE* outFile, Sequence* seq) {
  Mapping* best = bestMapping(0);
  if (best != NULL) {
    allBasesToLetters();

    uint32_t seqIndex = seq->sequenceIndex(best->start());
    uint32_t start = seq->positionInSeq(best->start());

    fprintf(outFile, "%s\t%d\t%s\t%d\t%d\t%s\t%c\t%d\t%d\t%s\t%s\n", id_,
            best->isComplement() ? 16 : 0, seq->info(seqIndex), start + 1, (uint32_t)best->score(),
            best->cigar(), '*', 0, 0, data_, quality_);
  }
}
