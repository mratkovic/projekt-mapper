/*
 * SuffixArray.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SUFFIXARRAY_H_
#define SUFFIXARRAY_H_

#include "Gene.h"
#include <vector>
#include <utility>

typedef unsigned long ulint;

class SuffixArray {

private:
	ulint _size;
	ulint* _array;
	char* _sequence;
	int _symbolNum;

	bool construct();
	void constructArray(ulint* str, ulint* SA, ulint n, int numOfSymbols);
	void radixSortPass(ulint* in, ulint* out, ulint* sequence, ulint n, int numOfSymbols);
	int compare(char *pattern, size_t len, ulint startPosition);


public:

	SuffixArray() {
		_sequence = 0;
		_array = 0;
		_size = 0;
		_symbolNum = 0;
	}

	SuffixArray(char* sequence, size_t size, int symbols) {
		_sequence = sequence;
		_size = size;
		_array = (ulint *) malloc(_size * sizeof(ulint));
		_symbolNum = symbols;
		construct();
	}

	~SuffixArray() {
		clear();
	}

	void clear() {
		if (_array) {
			free(_array);
		}
		if(_sequence){
			free(_sequence);
		}
		_size = 0;
	}

	ulint findStartingPositions(char *pattern, size_t patternLen, ulint id, std::vector<std::pair<int, ulint> > d);
	void printSuffixArray(FILE *out);
	void constructFromGene(bioutil::Gene* gene);

	ulint getSize() const {
		return _size;
	}
};
#endif /* SUFFIXARRAY_H_ */
