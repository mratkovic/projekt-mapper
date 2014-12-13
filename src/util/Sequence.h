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

	const char* data() const {
		return _data;
	}
	const char data(size_t i) const {
		return _data[i];
	}
	const char* description() const {
		return _description;
	}

	const int dataSize() const {
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
	size_t printSequence(FILE* outputFilePointer, int width = 80);
	bool readSequenceFromFASTA(FILE* inputFilePointer);
	void turnBaseToInt(char* array);
};
} /* namespace bioutil */

#endif /* SEQUENCE_H_ */
