/*
 * Mapping.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#include <cstdlib>
#include <cstdio>
#include "Mapping.h"

namespace bioutil {

Mapping::Mapping() {
	this->_score = this->_positionStart = this->_positionEnd = 0;
	this->_isComplement = false;
}
Mapping::Mapping(double score, size_t posStart, size_t posEnd, bool complemented = false) {
	this->_score = score;
	this->_positionStart = posStart;
	this->_positionEnd = posEnd;
	this->_isComplement = complemented;
}
char* Mapping::print() {
	char str[150];
	sprintf(str, "%5f, %lu, %lu, %d", this->_score, this->_positionStart, this->_positionEnd,
			this->_isComplement ? 1 : 0);
	return str;
}
Mapping::~Mapping() {
}

} /* namespace bioutil */
