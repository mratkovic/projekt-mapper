/*
 * SuffixArray.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SUFFIXARRAY_H_
#define SUFFIXARRAY_H_

#include <vector>
#include <utility>
#include "Sequence.h"


class SuffixArray {

private:
	int _size;
	int* _array;
	char* _sequence;
	int _symbolNum;

	bool construct();
	void constructArray(int* str, int* SA, int n, int numOfSymbols);
	void radixSortPass(int* in, int* out, int* sequence, int n, int numOfSymbols);
	int compare(char *pattern, size_t len, int startPosition);


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
		_array = (int *) malloc(_size * sizeof(int));
		_symbolNum = symbols;
		construct();
	}
	SuffixArray(FILE *file);

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

	int findStartingPositions(char *pattern, size_t patternLen, int id, std::vector<std::pair<int, int> > &d);
	void printSuffixArray(FILE *out);
	void constructFromSequence(bioutil::Sequence* gene);

	int getSize() const {
		return _size;
	}
};
#endif /* SUFFIXARRAY_H_ */
