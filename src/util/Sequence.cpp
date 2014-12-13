/*
 * Sequence.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include "Sequence.h"

#include "../util/BioUtil.h"
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>

#include <unistd.h>

#include "../util/UtilityFunctions.h"

namespace bioutil {

bool Sequence::readSequenceFromFASTA(FILE* inputFilePointer) {
	static char buffer[1024];
	this->clear();

	size_t dataLen = 0;
	size_t descriptionLen = 0;

	size_t numOfGetsForDescription = 0;
	size_t numOfGetsForData = 0;

	fpos_t prevPos;
	fgetpos(inputFilePointer, &prevPos);

	bool done = false;

	while (!done) {
		if (!fgets(buffer, sizeof buffer, inputFilePointer)) {
			assert(dataLen == 0 && descriptionLen == 0);
			return false;
		} else if (descriptionLen == 0) {
			assert(buffer[0] == '>');
		}

		int len = strlen(buffer);
		if (buffer[len - 1] == '\n') {
			done = true;
		}
		len = trimEnd(buffer);
		descriptionLen += len;
		++numOfGetsForDescription;
	}

	done = false;
	while (!done) {
		if (!fgets(buffer, sizeof buffer, inputFilePointer)) {
			break;
		} else if (buffer[0] == '>') {
			break;
		}
		dataLen += trimEnd(buffer);

		++numOfGetsForData;
	}

	this->_description = (char*) malloc(descriptionLen + 1);
	assert(this->_description);

	this->_data = (char*) malloc(dataLen + 1);
	assert(this->_data);
	this->_dataLen = dataLen;

	size_t positionDescription = 0;
	fsetpos(inputFilePointer, &prevPos);
	for (size_t i = 0; i < numOfGetsForDescription; ++i) {
		fgets(buffer, sizeof buffer, inputFilePointer);
		int len = trimEnd(buffer);
		memcpy(this->_description + positionDescription, buffer, len);
		positionDescription += len;
	}

	size_t positionData = 0;
	for (size_t i = 0; i < numOfGetsForData; ++i) {
		fgets(buffer, sizeof buffer, inputFilePointer);
		int len = trimEnd(buffer);

		memcpy(this->_data + positionData, buffer, len);
		positionData += len;
	}

	this->_data[dataLen] = 0;
	this->_description[descriptionLen] = 0;

	assert(positionData == dataLen);
	assert(positionDescription == descriptionLen);

	for (size_t i = 0; i < dataLen; ++i) {
		this->_data[i] = toupper(this->_data[i]);
		assert(isValidBaseACGT(this->_data[i]));
	}

	return dataLen > 0;
}

size_t Sequence::printSequence(FILE* outputFilePointer, int width) {
	size_t printed = 0;

	printed += fprintf(outputFilePointer, "%s\n", this->_description);

	for (size_t i = 0; i < this->dataSize(); i += width) {
		for (size_t j = i; j < std::min<int>(this->dataSize(), i + width); ++j) {
			printed += fprintf(outputFilePointer, "%c", this->data(j));
		}
		printed += fprintf(outputFilePointer, "\n");
	}

	return printed;
}

void Sequence::turnBaseToInt(char* array) {
	for (int i = 0; i < this->_dataLen; ++i) {
		array[i] = baseToInt(this->_data[i]);
	}
}


} /* namespace bioutil */

