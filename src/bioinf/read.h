/*
 * read.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SRC_CORE_READ_H_
#define SRC_CORE_READ_H_

#include <zlib.h>
#include <set>

#include "bioinf/sequence.h"
#include "external/kseq.h"
#include "core/suffix_array.h"
#include "bioinf/position.h"
#include "util/utility_functions.h"

#define MAX_KEEP 80
#define KEEP_FACTOR 1.2f

KSEQ_INIT(int, read)

namespace bioinf {

class Read {
 public:
  Read(float keepRatio = KEEP_FACTOR, uint32_t maxPositions = MAX_KEEP);
  ~Read();

  void clear();
  bool readNextFromFASTQ(kseq_t *seq);
  void allBasesToSmallInt();
  void allBasesToLetters();

  Read* getReverseComplement();
  void printReadSAM(FILE* outFile, Sequence* seq);

  const char* data();
  const char* id();
  uint32_t dataLen();

  void addPosition(uint32_t score, uint32_t start, uint32_t end,
                   bool isComplement = false, const char* cigar = NULL,
                   uint32_t cigarLen = 0, uint32_t secondaryScore = 0);
  void addPosition(Position* position);

  std::set<Position*, ptr_compare<Position> >& positions();
  Position* bestPosition(uint32_t index);
  uint32_t positionsSize();
  void keepRatio(float keepRatio);
  void maxPositions(uint32_t maxPositions);
  bool basesInt();

 private:
  char* id_;
  char* data_;
  uint32_t dataLen_;

  char* optional_identifier_;
  char* quality_;

  bool basesAsInt_;
  std::set<Position*, ptr_compare<Position> > positions_;

  float keepRatio_;
  uint32_t maxPositions_;

};

}

#endif /* SRC_CORE_READ_H_ */
