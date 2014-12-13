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
#define KEEP_FACTOR 2

namespace bioutil {

Read::Read() {
	_id = 0;
	_data = 0;
	_dataLen = 0;
}
Read::~Read() {
	printf("Dest\n");
	this->clean();
}

void Read::addMapping(Mapping m) {
	this->_mappings.insert(m);
	// update set
	while (this->_mappings.size() > 0) {
		double bestScore = this->_mappings.rbegin()->getScore();
		double worstScore = this->_mappings.begin()->getScore();

		double ratio = bestScore / worstScore;
		if (ratio > KEEP_FACTOR) {
			this->_mappings.erase(this->_mappings.begin());
		} else {
			break;
		}
	}
}
Mapping Read::getBestMapping() {
	if (this->_mappings.size() > 0) {
		return *this->_mappings.rbegin();
	}
	Mapping m(-1, 0, 0, false);
	return m;

}
std::multiset<Mapping> Read::getMappings() {
	return this->_mappings;
}

void Read::getAllReadsFromFASTQ(FILE* input, std::vector<Read*> &reads) {
	Read tmpRead;
	while (tmpRead.readFromFASTQ(input)) {

		Read* read = new Read();
		*read = tmpRead;
		reads.push_back(read);
	}
}

void Read::clean() {
	if (this->_data) {
		free(this->_data);
	}
	if (this->_id) {
		free(this->_id);
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
	if (this->_mappings.size() > 0)
		printed += fprintf(outputFilePointer, "%s\; %s", this->_id, this->getBestMapping().print());

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
