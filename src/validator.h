/*
 * validator.h
 *
 *  Created on: Jan 4, 2015
 *      Author: marko
 */

#ifndef SRC_VALIDATOR_H_
#define SRC_VALIDATOR_H_

#include "read.h"
#include <cstdio>

namespace bioutil {

class Validator {
public:
	static void validateSAM(FILE* ref, FILE* test);
	static void validateWGSIM(Read* read);
};

} /* namespace bioutil */

#endif /* SRC_VALIDATOR_H_ */
