/*
 * validator.cpp
 *
 *  Created on: Jan 4, 2015
 *      Author: marko
 */

#include "validator.h"
#include <map>
#include <cstdio>
#include <stdint.h>

namespace bioutil {

void Validator::validateSAM(FILE* ref, FILE* test) {
	std::map<char*, uint32_t> reads;
	// TODO
}
void Validator::validateWGSIM(Read* read) {
	//wgsim header
	//readID_posStart_posEnd_.....
	// TODO

}

} /* namespace bioutil */
