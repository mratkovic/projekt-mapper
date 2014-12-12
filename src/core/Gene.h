/*
 * Gene.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef GENE_H_
#define GENE_H_

#include <cstring>
#include <string>
#include <vector>
#include <cstdlib>

namespace bioutil {

typedef unsigned long ulint;

class Gene {
private:
	char* _data;
	ulint _dataLen;
	char* _description;

public:
	Gene() {
		_description = _data = 0;
		_dataLen = 0;
	}
	~Gene() {
		clear();
	}

	const char* data() const {
		return _data;
	}
	const char data(size_t i) const {
		return _data[i];
	}
	const char* description() const {
		return _description;
	}

	const ulint dataSize() const {
		return _dataLen;
	}

	void clear() {
		if (_data) {
			free(_data);
		}
		if (_description) {
			free(_description);
		}
		_dataLen = 0;
	}
	size_t printGene(FILE* outputFilePointer, int width = 80);
	bool readGeneFromFASTA(FILE* inputFilePointer);
	void turnBaseToInt(char* array);
};
} /* namespace bioutil */

#endif /* GENE_H_ */
