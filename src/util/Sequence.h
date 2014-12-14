/*
 * Sequence.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <cstring>
#include <string>
#include <vector>
#include <cstdlib>

namespace bioutil {

class Sequence {
private:
	char* _data;
	int _dataLen;
	char* _description;

public:
	Sequence() {
		_description = _data = 0;
		_dataLen = 0;
	}
	~Sequence() {
		clear();
	}

	char* data() {
		return _data;
	}
	char data(size_t i) {
		return _data[i];
	}
	char* description() {
		return _description;
	}

	int dataSize() {
		return _dataLen;
	}

	void clear() {
		if (_data) {
			free(_data);
			_data = 0;
		}
		if (_description) {
			free(_description);
			_description = 0;
		}
		_dataLen = 0;
	}
	size_t printSequence(FILE* outputFilePointer, int width = 80);
	bool readSequenceFromFASTA(FILE* inputFilePointer);
	void turnBaseToInt();
};
} /* namespace bioutil */

#endif /* SEQUENCE_H_ */
