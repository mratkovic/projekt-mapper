/*
 * mapper.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#ifndef SRC_UTIL_MAPPER_H_
#define SRC_UTIL_MAPPER_H_

#include <stdint.h>
#include <omp.h>

#include "suffix_array.h"

#include "../bioutil/read.h"
#include "../bioutil/sequence.h"
#include "../util/fenwick.h"

#define KMER_K 10
#define WINDOW_SIZE 2

#define SW_START_OFFSET 20
#define SW_END_OFFSET 20

using namespace bioutil;

class Mapper {

 private:
  static void runLIS(int startIndex, int endIndex,
                     std::vector<std::pair<uint32_t, uint32_t> > &pos,
                     Read* read);
  static void runLCSk(int startIndex, int endIndex,
                      std::vector<std::pair<uint32_t, uint32_t> > &pos,
                      Read* read);

  static void mapReadToSuffixArray(Read* read, SuffixArray* sa,
                                   bool generateCIGAR);
  static void getKmerPositions(
      Read* read, SuffixArray* sa,
      std::vector<std::pair<uint32_t, uint32_t> > &positions, int kmerStart);
  static void fillMappings(Read* read, SuffixArray* sa);
  static void getPositions(Read* read, SuffixArray* sa,
                           bool complement);

 public:
  static void mapAllReads(char* readsInPath, char* solutionOutPath,
                          SuffixArray* sa, Sequence* seq,
                          bool generateCIGAR);

};
#endif
