/*
 * Read.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SRC_CORE_READ_H_
#define SRC_CORE_READ_H_

#include <cstdlib>
#include <cstdio>
#include <vector>

namespace bioutil {

class Read {
private:
	char* _id;
	char* _data;
public:
	Read();
	virtual ~Read();
	bool readFromFASTQ(FILE *inputFile);

};

static void getAllReadsFromFASTQ(FILE * input, std::vector <Read*> reads);
}
/* namespace bioutil */

#endif /* SRC_CORE_READ_H_ */
