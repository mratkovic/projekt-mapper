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

#define KMER_K 20
#define WINDOW_SIZE 2

// TODO trenutno bzvz stoji samo
#define MAX_NUMBER_OF_KMER_POSITIONS 1000000

#define SSW_MATCH 3
#define SSW_MISMATCH 1
#define SSW_GAP_OPEN 3
#define SSW_GAP_EXTEND 2

#define MAX_EDIT_DIST_FACTOR 0.1l

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
