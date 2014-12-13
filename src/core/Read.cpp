/*
 * Read.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include <algorithm>
#include <cassert>
#include <cstring>

#include "Read.h"
#include "../util/UtilityFunctions.h"

#define LINE_SIZE 30000

namespace bioutil {

Read::Read() {
	_id = 0;
	_data = 0;
	_mappingInfo = 0;
}
Read::~Read() {
	if (_id) {
		free(_id);
	}
	if (_data) {
		free(_data);
	}
	if (_mappingInfo) {
		free(_mappingInfo);
	}
}

void Read::setMappingInfo(char* info) {
	size_t len = strlen(info) + 1;
	this->_mappingInfo = (char *) malloc(len * sizeof(char));
	strcpy(this->_mappingInfo, info);
}
void Read::getAllReadsFromFASTQ(FILE* input, std::vector<Read*> &reads) {
	Read tmpRead;
	printf("Ucitavam readove\n");
	while (tmpRead.readFromFASTQ(input)) {

		Read* read = new Read();
		*read = tmpRead;
		reads.push_back(read);
	}
}
bool Read::readFromFASTQ(FILE* inputFile) {
	char id[LINE_SIZE + 1];
	char data[LINE_SIZE + 1];
	char plus[LINE_SIZE + 1];
	char quality[LINE_SIZE + 1];
	id[LINE_SIZE] = data[LINE_SIZE] = plus[LINE_SIZE] = quality[LINE_SIZE] = 0;

	if (!fgets(id, sizeof id, inputFile)) {
		return 0;
	}
	assert(id[LINE_SIZE] == 0);

	if (!fgets(data, sizeof data, inputFile)) {
		return 0;
	}
	assert(data[LINE_SIZE] == 0);

	assert(fgets(plus, sizeof plus, inputFile));
	assert(plus[LINE_SIZE] == 0);
	assert(fgets(quality, sizeof quality, inputFile));
	assert(quality[LINE_SIZE] == 0);

	int lenID = trim(id);
	this->_id = (char *) malloc((lenID + 1) * sizeof(char));
	for (int i = 0; i < lenID; ++i) {
		this->_id[i] = id[i];
	}
	this->_id[lenID] = 0;

	int lenData = trim(data);
	this->_data = (char *) malloc((lenData + 1) * sizeof(char));
	for (int i = 0; i < lenData; ++i) {
		this->_data[i] = data[i];
	}
	this->_data[lenData] = 0;
	this->_dataLen = lenData;
	return lenID > 0 && lenData > 0;
}

size_t Read::printRead(FILE* outputFilePointer, int width) {
	size_t printed = 0;
	if (this->_mappingInfo)
		printed += fprintf(outputFilePointer, "%s\n", this->_mappingInfo);

	for (size_t i = 0; i < this->getDataLen(); i += width) {
		for (size_t j = i; j < std::min<int>(this->getDataLen(), i + width); ++j) {
			printed += fprintf(outputFilePointer, "%c", this->getData(j));
		}
		printed += fprintf(outputFilePointer, "\n");
	}

	return printed;
}
}
/* namespace bioutil */
