/*
 * mapping.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#ifndef SRC_UTIL_MAPPING_H_
#define SRC_UTIL_MAPPING_H_

#include <stdint.h>

namespace bioutil {

class Mapping {
private:
	double score_;
	uint32_t start_;
	uint32_t end_;
	bool complement_;
	char* cigar_;

public:
	Mapping();
	Mapping(double score, uint32_t start, uint32_t end, bool complemented, const char* cigar, uint32_t cigarLen);
	virtual ~Mapping();

	bool isComplement();

	double score() const;

	uint32_t end();
	uint32_t start();

	void end(uint32_t end);
	void start(uint32_t start);

	const char* cigar();
	void cigar(const char* cigar, uint32_t len);
	void setComplement(bool complemented);
	bool operator <(const Mapping &other) const{
			return score_ <= other.score();
		}

};

} /* namespace bioutil */

#endif /* SRC_UTIL_MAPPING_H_ */
