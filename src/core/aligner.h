/*
 * aligner.h
 *
 *  Created on: May 9, 2015
 *      Author: marko
 */

#ifndef SRC_CORE_ALIGNER_H_
#define SRC_CORE_ALIGNER_H_

#include <core/solver.h>
#include <core/lcsk_solver.h>
#include <core/incremental_lcsk_solver.h>

#include <algorithm>
#include <omp.h>

#define THREAD_NUM std::max(omp_get_num_procs() / 2, 12)

class Aligner {
 public:
  Aligner(uint32_t kmerK = KMER_K, uint32_t lowerMatchesLimit = MIN_MATCH_NUM,
          uint32_t upperMatchesLimit = MAX_MATCH_NUM, uint32_t threadNum =
              THREAD_NUM, float keepFactor = KEEP_FACTOR, uint32_t maxPos =MAX_KEEP);

  virtual ~Aligner();
  void run(int argc, char** argv);
  void constructSA(char* fastaInPath, char* saOutputPath);
  void mapReads(char* fastaInPath, char* saInputPath, char* readsInPath, char* SAMOutputPath);

 private:

  void mapReads(int argc, char** argv);
  void verboseUsageAndExit();

  uint32_t kmerK_;
  uint32_t lowerMatchesLimit_;
  uint32_t upperMatchesLimit_;
  uint32_t threadNum_;
  float keepFactor_;
  uint32_t maxPos_;
};

#endif /* SRC_CORE_ALIGNER_H_ */
