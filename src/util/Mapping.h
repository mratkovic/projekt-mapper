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
	size_t _positionStart;
	size_t _positionEnd;
	bool _isComplement;

public:
	Mapping();
	Mapping(double score, size_t posStart, size_t posEnd, bool complemented);
	virtual ~Mapping();
	char* print();
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

	size_t getPositionEnd() const {
		return _positionEnd;
	}

	size_t getPositionStart() const {
		return _positionStart;
	}
};

} /* namespace bioutil */

#endif /* SRC_UTIL_MAPPING_H_ */
