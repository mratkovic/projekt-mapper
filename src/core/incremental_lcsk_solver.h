/*
 * incremental_lcsk_solver.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_CORE_INCREMENTAL_LCSK_SOLVER_H_
#define SRC_CORE_INCREMENTAL_LCSK_SOLVER_H_

#define MAX_MATCH_NUM 45
#define MIN_MATCH_NUM 20

#include <cstdint>
#include <vector>
#include "util/util_structures.h"
#include "bioinf/sequence.h"
#include "lcsk_solver.h"
#include "suffix_array.h"

class IncrementalLCSkSolver : public LCSkSolver {

 public:
  IncrementalLCSkSolver(bioinf::Sequence* seq);
  virtual ~IncrementalLCSkSolver();
  virtual void printInfo();

  uint32_t maxMatchNum_;
  uint32_t minMatchNum_;

 private:
  virtual void runLCSkpp(int startIndex, int endIndex,
                         std::vector<triplet_t<uint32_t>> &pos,
                         bioinf::Read* read);
  virtual void fillPositions(bioinf::Read* read);
  virtual uint32_t getKmerPositions(bioinf::Read* read,
                                    std::vector<triplet_t<uint32_t>> &positions,
                                    int kmerStart, uint32_t inittialLen);

};

#endif /* SRC_CORE_INCREMENTAL_LCSK_SOLVER_H_ */
