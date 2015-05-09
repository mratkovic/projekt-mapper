/*
 * lcsk_solver.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_CORE_LCSK_SOLVER_H_
#define SRC_CORE_LCSK_SOLVER_H_

#include <cstdint>
#include <vector>
#include "util/util_structures.h"
#include "bioinf/sequence.h"
#include "solver.h"
#include "suffix_array.h"

#define KMER_K 15
#define WINDOW_SIZE 2

#define SSW_MATCH 5
#define SSW_MISMATCH 4
#define SSW_GAP_OPEN 8
#define SSW_GAP_EXTEND 8


class LCSkSolver : public Solver {
 public:
  LCSkSolver(bioinf::Sequence* seq);
  virtual ~LCSkSolver();

  virtual void findReadPosition(bioinf::Read* read);
  void readSuffixArrayFromFile(const char* saInPath);
  virtual void printInfo();


  bioinf::Sequence* seq();
  SuffixArray* sa();

  uint32_t kmerK_;
  uint32_t windowSize_;

 private:
  virtual void runLCSkpp(int startIndex, int endIndex,
                 std::vector<std::pair<uint32_t, uint32_t>> &pos,
                 bioinf::Read* read);
  virtual void fillPositions(bioinf::Read* read);
  virtual void getKmerPositions(bioinf::Read* read,
                        std::vector<std::pair<uint32_t, uint32_t>> &positions,
                        int kmerStart);

  bioinf::Sequence* seq_;
  SuffixArray* sa_;
};

#endif /* SRC_CORE_LCSK_SOLVER_H_ */
