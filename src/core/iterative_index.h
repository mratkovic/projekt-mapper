/*
 * iterative_index.h
 *
 *  Created on: Feb 16, 2016
 *      Author: marko
 */

#ifndef SRC_CORE_ITERATIVE_INDEX_H_
#define SRC_CORE_ITERATIVE_INDEX_H_

#include <core/index.h>

class IterativeIndex : public Index {
 public:
    IterativeIndex(){};
    virtual ~IterativeIndex(){};

    virtual const int* search(const char* pattern, int length, int* numOfSolutions)=0;

    virtual const int* iterativeSearch(const char* pattern, int length, int startLen,
                               int* numOfSolutions, int* finalLen, const int solUpperLimit,
                               const int solLowerLimit)=0;

    virtual void save(FILE *out)=0;
};


#endif /* SRC_CORE_ITERATIVE_INDEX_H_ */
