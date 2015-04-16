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

#include "bioinf/read.h"
#include "bioinf/sequence.h"
#include "solver.h"

#define MAX_TMP_NAME_LEN 50

using namespace bioinf;

class Mapper {

 public:

  Mapper(Sequence* seq, Solver* solver, uint32_t threadNum = omp_get_num_procs(), float minKeepScoreRatio = KEEP_FACTOR, uint32_t maxPositionsPerRead = MAX_KEEP);
  void mapAllReads(char* readsInPath, char* solutionOutPath);

 private:
  // TODO REMOVE OR REFACTOR
  static void runLIS(int startIndex, int endIndex,
                     std::vector<std::pair<uint32_t, uint32_t> > &pos,
                     Read* read);

  void fillSAMHeader(FILE* out);

  void copyFromTmpFileAndDelete(char* tmpFileName, FILE* src, FILE *dest);
  void createTmpFiles(FILE* tempFiles[], char tmpFileNames[][MAX_TMP_NAME_LEN],
                      uint32_t numOfFiles);
  void mergeTmpFiles(char fileNames[][MAX_TMP_NAME_LEN], FILE* tmpFiles[],
                     FILE* solutionFile, int numberOfFiles);

  Sequence* seq_;
  Solver* solver_;
  uint32_t threadNum_;

  float minKeepScoreRatio_;
  uint32_t maxPositionsPerRead_;
};
#endif
