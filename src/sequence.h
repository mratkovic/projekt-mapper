/*
 * sequence.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_


#include <cstdlib>
#include <cstdio>
#include <stdint.h>

namespace bioutil {

class Sequence {

private:
	char* data_;
	uint32_t dataLen_;
	char* info_;
	char* comment_;

public:
	Sequence();
	Sequence(char* data, uint32_t dataLen, char* info, char* comment);
	~Sequence();
	void clear();
	void readSequenceFromFASTA(FILE* fastaIn);
	const char* data();
	uint32_t dataLen();
	const char* info();
	const char* comment();

};
}
#endif /* SEQUENCE_H_ */
