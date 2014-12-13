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
#include <set>
#include "SuffixArray.h"
#include "Mapping.h".h"
namespace bioutil {

class Read {
private:
	char* _id;
	char* _data;
	size_t _dataLen;
	std::multiset<Mapping> _mappings;

public:
	Read();
	virtual ~Read();
	void clean();
	bool readFromFASTQ(FILE *inputFile);
	static void getAllReadsFromFASTQ(FILE * input, std::vector<Read*> &reads);
	size_t printRead(FILE* outputFilePointer, int width = 80);

	char getData(int i) const {
		return _data[i];
	}
	char* data() const {
		return _data;
	}

	char* getId() const {
		return _id;
	}

	size_t getDataLen() const {
		return _dataLen;
	}

	void addMapping(Mapping m);
	Mapping getBestMapping();
	std::multiset<Mapping> getMappings();

};

}
/* namespace bioutil */

#endif /* SRC_CORE_READ_H_ */
