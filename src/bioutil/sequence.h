/*
 * sequence.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <vector>

#include <cstdlib>
#include <cstdio>
#include <cstdint>

namespace bioutil {

class Sequence {
 public:
  Sequence();
  //Sequence(char* data, uint32_t dataLen, char* info, char* comment);
  ~Sequence();
  void clear();
  void readSequencesFromFASTA(FILE* fastaIn);
  void readSingleSequenceFromFASTA(FILE* fastaIn);
  void allBasesToSmallInt();
  void allBasesToLetters();
  const char* data();
  uint32_t dataLen();
  const char* info(uint32_t index);
  uint32_t positionInSeq(uint32_t positionGlobal);
  uint32_t sequenceIndex(uint32_t positionGlobal);
  uint32_t numOfSequences();
  uint32_t seqLen(uint32_t index);

  bool basesInt();

 private:
  char* data_;
  uint32_t dataLen_;

  uint32_t numOfSequences_;
  std::vector<uint32_t> seqEndIndex_;
  std::vector<char*> info_;

  bool basesAsInt;

};
}
#endif /* SEQUENCE_H_ */
