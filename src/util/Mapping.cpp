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
//	this->_editDistance = new int;
//	*this->_editDistance = -1;
}
Mapping::Mapping(double score, int posStart, int posEnd, bool complemented = false) {
	this->_score = score;
	this->_positionStart = posStart;
	this->_positionEnd = posEnd;
	this->_isComplement = complemented;
//	this->_editDistance = new int;
//	*this->_editDistance = -1;
}

void Mapping::fillDetails(char* printBuffer) {
	sprintf(printBuffer, "%f, %d, %d, %d", this->_score, this->_positionStart, this->_positionEnd,
			this->_isComplement ? 1 : 0);
}
Mapping::~Mapping() {
//	if (_editDistance) {
//		delete _editDistance;
//	}
}

} /* namespace bioutil */
