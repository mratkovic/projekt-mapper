/*
 * index.h
 *
 *  Created on: Feb 16, 2016
 *      Author: marko
 */

#ifndef SRC_CORE_INDEX_H_
#define SRC_CORE_INDEX_H_


class Index {
 public:
    Index(){};
    virtual ~Index(){};

    virtual const int* search(const char* pattern, int length, int* numOfSolutions)=0;
    virtual void save(FILE *out)=0;
};

#endif /* SRC_CORE_INDEX_H_ */
