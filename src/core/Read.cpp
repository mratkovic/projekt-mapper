/*
 * Read.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include <cassert>
#include <cstring>

#include "Read.h"
#include "../util/Util.h"
#define LINE_SIZE 30000

namespace bioutil {

Read::Read() {
	_id = 0;
	_data = 0;
}
Read::~Read() {
	if (_id) {
		free(_id);
	}
	if (_data) {
		free(_data);
	}
}
void getAllReadsFromFASTQ(FILE* input, std::vector<Read*> &reads) {
	Read tmpRead;
	while(tmpRead.readFromFASTQ(input)) {
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
	assert(data[100000] == 0);

	assert(fgets(plus, sizeof plus, inputFile));
	assert(plus[100000] == 0);
	assert(fgets(quality, sizeof quality, inputFile));
	assert(quality[100000] == 0);

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

	return lenID > 0 && lenData > 0;
}

}
/* namespace bioutil */
