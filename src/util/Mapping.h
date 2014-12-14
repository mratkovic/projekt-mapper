/*
 * Mapping.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#ifndef SRC_UTIL_MAPPING_H_
#define SRC_UTIL_MAPPING_H_
#include <cstdlib>

namespace bioutil {

class Mapping {
private:
	double _score;
	int _positionStart;
	int _positionEnd;
	bool _isComplement;

	//int* _editDistance;

public:
	Mapping();
	Mapping(double score, int posStart, int posEnd, bool complemented);
	virtual ~Mapping();
	void fillDetails(char* printBuffer);
	bool operator <(const Mapping &other) const{
		if (this->_score < other.getScore()) {
			return true;
		}
		return false;
	}

	bool isComplement() const {
		return _isComplement;
	}

	double getScore() const {
		return _score;
	}

	int getPositionEnd() const {
		return _positionEnd;
	}

	int getPositionStart() const {
		return _positionStart;
	}
};

} /* namespace bioutil */

#endif /* SRC_UTIL_MAPPING_H_ */
