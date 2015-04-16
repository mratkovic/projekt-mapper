/*
 * read.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include "read.h"
#include "util/utility_functions.h"
#include "bioutil.h"

#define LINE_SIZE 30000

namespace bioinf {

Read::Read(float keepRatio, uint32_t maxPositions)
    : keepRatio_(keepRatio),
      maxPositions_(maxPositions) {
  id_ = data_ = optional_identifier_ = quality_ = 0;
  dataLen_ = 0;
  basesAsInt_ = false;

}
Read::~Read() {
  clear();
}

void Read::addPosition(uint32_t score, uint32_t start, uint32_t end,
                       bool isComplement, const char* cigar,
                       uint32_t cigarLen, uint32_t secondaryScore) {

  Position* p = new Position(score,secondaryScore, start, end, isComplement, cigar, cigarLen);
  addPosition(p);

}

void Read::addPosition(Position* p) {
  positions_.insert(p);
  while (positions_.size() > 1
      && ((((double) (*positions_.rbegin())->score()
          / (*positions_.begin())->score()) > KEEP_FACTOR)
          || positions_.size() > MAX_KEEP)) {

    delete *positions_.begin();
    positions_.erase(positions_.begin());
  }
}

Position* Read::bestPosition(uint32_t index) {
  std::multiset<Position*, ptr_compare<Position> >::reverse_iterator it =
      positions_.rbegin();
  uint32_t cntr = 0;

  for (; it != positions_.rend(); ++it, ++cntr) {
    if (cntr == index) {
      return *it;
    }
  }
  return NULL;
}
std::set<Position*, ptr_compare<Position> >& Read::positions() {
  return positions_;
}
uint32_t Read::positionsSize() {
  return positions_.size();
}

const char* Read::data() {
  return data_;
}
const char* Read::id() {
  return id_;
}

void Read::keepRatio(float keepRatio) {
  keepRatio_ = keepRatio;
}
void Read::maxPositions(uint32_t maxPositions) {
  maxPositions_ = maxPositions;
}

uint32_t Read::dataLen() {
  return dataLen_;
}

bool Read::basesInt() {
  return basesAsInt_;
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

  uint32_t cntr = 0;
  for (auto it = positions_.begin(); it != positions_.end(); ++it, ++cntr) {
    delete *it;
  }
  dataLen_ = 0;

}

bool Read::readNextFromFASTQ(kseq_t* seq) {
  if (kseq_read(seq) < 0) {
    kseq_destroy(seq);
    return false;
  }

  id_ = seq->name.s;
  data_ = seq->seq.s;
  dataLen_ = seq->seq.l;
  quality_ = seq->qual.s;

  seq->name.s = NULL;
  seq->name.m = 0;
  seq->seq.s = NULL;
  seq->seq.m = 0;
  seq->qual.s = NULL;
  seq->qual.m = 0;

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
  Position* best = bestPosition(0);

  if (best != NULL) {

    allBasesToLetters();

    uint32_t seqIndex = seq->sequenceIndex(best->start());
    uint32_t start = seq->positionInSeq(best->start());

    fprintf(outFile, "%s\t%d\t%s\t%d\t%d\t%s\t%c\t%d\t%d\t%s\t%s\n", id_,
            best->isComplement() ? 16 : 0, seq->info(seqIndex), start + 1,
            (uint32_t) best->score(), best->cigar(), '*', 0, 0, data_,
            quality_);

//    if (positions_.size() > 4 && false) {
//      fprintf(stdout, "%s\t%d\t%s\t%d\t%d\n", id_,
//              best->isComplement() ? 16 : 0, seq->info(seqIndex), start + 1,
//              (uint32_t) best->score());
//    }

  } else {
    // TODO: nije mapiran
  }

}

}  // end namespace
