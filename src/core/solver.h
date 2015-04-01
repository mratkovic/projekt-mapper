/*
 * solver.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_CORE_SOLVER_H_
#define SRC_CORE_SOLVER_H_

#include "bioinf/read.h"

class Solver {
 public:
  Solver() {
  }
  virtual ~Solver() {
  }

  virtual void printInfo() = 0;
  virtual void findReadPosition(bioinf::Read* read) = 0;
};

#endif /* SRC_CORE_SOLVER_H_ */
